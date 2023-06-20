function [ fval, grad, H ] = corr_grad_hess(x, Nm, Nk, wmk, Lmi, Omki, mu_m, nu_k, Et)
%FP_UPDATE Summary of this function goes here

%   Detailed explanation goes here
ti = reshape(x(1:3), 3, 1);
ui = reshape(x(4:6), 3, 1);

DIMS = 3;
% Derive quaternion value
v = sqrt(ui' * ui);
q = [cos(v/2); 0.5*sinc( v/(2*pi) ) * ui];

% Generate rotation matrix
Ri = quat2rotm(q');

% Compute R*iOm
ROmi = mtimesx(Ri, Omki);

% Compute Smk
ROmiR = mtimesx(ROmi, Ri, 't'); % R * Omk^-1 * R'
Smki = reshape(repmat(Lmi, 1, Nk, 1), DIMS, DIMS, Nm*Nk) + repmat(ROmiR, 1, 1, Nm);
Smk = zeros(DIMS, DIMS,Nm*Nk);
det_Smk = zeros(Nm*Nk,1);
% det_Smki = zeros(Nm*Nk:,1);

for ix = 1:Nm*Nk
    Smk(:, :, ix) = inv(Smki(:, :, ix));
    det_Smk(ix) = det(Smk(:, :, ix));
%     det_Smki(ix) = det(Smki(:, :, ix));
end

% Compute ymk
Rnu = mtimesx(Ri, nu_k);
mu_Rnu = reshape(repmat(mu_m, Nk, 1, 1), 3, 1, Nm*Nk) - ...
         repmat(reshape(Rnu, 3, 1, Nk), 1, 1, Nm);
ymk = mu_Rnu - repmat(ti,1,1,Nm*Nk);

% Compute Symk
Symk = mtimesx(Smk, ymk);

PI_32 = 1/(2*pi)^(1.5);
% Compute fmk .* (2*pi)^(-3/2)
fmk = reshape(wmk, 1, 1, Nm*Nk)  .* PI_32 .* reshape(sqrt(det_Smk), 1,1, Nm*Nk)...
       .* exp(-0.5 * mtimesx(ymk, 't', Symk));
fmk = -fmk / Et;

% Evaluate correlation function
fval = sum(fmk, 3);

% Evaluate gradient
if nargout > 1

    fmk_vec = reshape(repmat(fmk, 3, 1, 1), 3, 1, Nm*Nk);
    fmk_mat = reshape(repmat(fmk, 9, 1, 1), 3, 3, Nm*Nk);

    % Compute partial derivative wrt t
    dFdt = sum(fmk_vec .* Symk, 3);

    % Compute partial derivative wrt R
    nu_mk = repmat(reshape(nu_k, 3, 1, Nk), 1, 1, Nm);
    ROmi_mk = repmat(ROmi, 1, 1, Nm);
    Synut_mk = mtimesx(Symk, nu_mk, 't');
    Sy_ytSmk = mtimesx(Symk, Symk, 't');
    Sy_ytSROmi_mk = mtimesx(Sy_ytSmk, ROmi_mk);
    SmkROmi = mtimesx(Smk, repmat(ROmi, 1, 1, Nm    ));

    dFdR = sum(fmk_mat .* (Synut_mk + Sy_ytSROmi_mk - SmkROmi), 3);

    grad(1:3) = dFdt;
    if nargout == 2
        dRdu = partial_wrt_u(ui);
    else
        [dRdu, d2Rdu2] = partial_wrt_u(ui);
    end
    for ix = 1:3
        grad(3+ix) = trace(dFdR' * dRdu(:, :, ix));
    end
end

% Evaluate Hessian
if nargout > 2

    % Hessian matrix
    H = zeros(6,6);

    % Compute dt^2 terms
    H(1:3, 1:3) = sum(fmk_mat .* (Sy_ytSmk - Smk), 3);

    Delta = zeros(3,3,3);
    Gamma = zeros(3,3,3);

    % Compute dtdr and drdt terms
    % trace indices
    tr_inds = bsxfun(@plus, (0:Nm*Nk-1)*9, uint32([1 5 9]'));

    OmiRtSmk = permute(SmkROmi, [2 1 3]);
    OmiRtSy_ytSmk = permute(Sy_ytSROmi_mk, [2 1 3]);

    for ix = 1:3
        ROmiAt = mtimesx(ROmi, dRdu(:, :, ix)');
        AOmiRt = permute(ROmiAt, [2 1 3]);

        Za = ROmiAt + AOmiRt;
        Za_mk = repmat(Za, 1, 1, Nm);

        Anu_k = mtimesx(dRdu(:, :, ix), nu_k);
        Anu_mk = repmat(reshape(Anu_k, 3, 1, Nk), 1, 1, Nm);

        ZaSymk = mtimesx(Za_mk, Symk);
        SmkROmiAt = mtimesx(Smk, repmat(ROmiAt, 1, 1, Nm));

        da_mk = -sum( SmkROmiAt(tr_inds) );
        da_mk = repmat(reshape(da_mk, 1, 1, Nm*Nk), 3, 1, 1);
        da_mat = repmat(da_mk, 1, 3, 1);

        qa = mtimesx(ZaSymk + 2 * Anu_mk, 't', Symk);
        qa_mk = repmat(qa, 3, 1, 1);
        qa_mat = repmat(qa, 3,3,1);

        S__ZaSy_Anu_mk = mtimesx(Smk, ZaSymk + Anu_mk);

        % dtdr
        H(1:3, 3+ix) = sum(fmk_vec .* ( (0.5*qa_mk + da_mk).*Symk ...
                                        - S__ZaSy_Anu_mk), 3);

        % dr^2
        OmAt = repmat(mtimesx(Omki, dRdu(:, :, ix)'), 1, 1, Nm);

        OmiRtSZa_mk = mtimesx(OmiRtSmk, Za_mk);
        %Dba = da_mat .* OmiRtSmk + mtimesx(-OmAt + OmiRtSZa_mk, Smk);
        Dba = mtimesx(-OmAt + OmiRtSZa_mk, Smk);

        da_qa = (2*da_mat + qa_mat);
        da_qa__db = -da_qa .* OmiRtSmk;
        nuytSmk = permute(Synut_mk, [2 1 3]);
        da_qa__qb = da_qa .* ( OmiRtSy_ytSmk + nuytSmk );
        ZaSyyt = mtimesx(ZaSymk, ymk, 't');
        yytSZa = permute(ZaSyyt, [2 1 3]);
        % T
        nunutAt_k = mtimesx(reshape(nu_k, 3, 1, Nk), ...
                            reshape(Anu_k, 3, 1, Nk), 't');
        nunutAt_mk = repmat(nunutAt_k, 1, 1, Nm);
        ytnutAt_mk = mtimesx(ymk, Anu_mk, 't' );
        Anyt_mk = permute(ytnutAt_mk, [2 1 3]);


        dqadrb_b =  -2 * mtimesx(nu_mk, mtimesx(ZaSymk, 't', Smk)) ...
                   - 2 * mtimesx(OmiRtSmk, mtimesx(ZaSyyt + yytSZa, Smk)) ...
                   + 2 * mtimesx(OmAt, Sy_ytSmk) ...
                   - 2 * mtimesx(OmiRtSmk, mtimesx(ytnutAt_mk + Anyt_mk, Smk))...
                   - 2 * mtimesx(nunutAt_mk, Smk);

        Delta_mk = 2*Dba + da_qa__db + da_qa__qb + dqadrb_b;
        Delta(:, :, ix) = 0.5 * sum(fmk_mat .* Delta_mk, 3);

        Gamma_mk = -2 * OmiRtSmk + 2 * OmiRtSy_ytSmk ...
                   +2 * nuytSmk;
        Gamma(:, :, ix) = 0.5 * sum(fmk_mat .* Gamma_mk, 3);
    end
    for ix = 1:3
        for kx = 1:3
            H(3+ix, 3+kx) = trace(Delta(:, :, kx) * dRdu(:, :, ix)) + ...
                            trace(Gamma(:, :, kx) * d2Rdu2(:, :, ix, kx));
        end
    end
%      trace(d2Rdu2(:, :, ix, kx) * dRdu(:, :, 1) + d2Rdu2(:, :, 1, 1)
    H(4:6, 1:3) = H(1:3, 4:6)';

end

end

function [x_opt, score, grad, Hess] = gmm_registration(gmm_source, gmm_target, x_init)

     % Decompose source gmm
     Ns = uint32(gmm_source.n_components);
     w_s = transpose(gmm_source.weights);
     S_s = reshape(transpose(gmm_source.covs), Ns, 3, 3);
     S_s = permute(S_s, [3, 2, 1]);
     mu_s = gmm_source.means;

     % Decompose target gmm
     Nt = uint32(gmm_target.n_components);
     w_t = transpose(gmm_target.weights);
     S_t = reshape(transpose(gmm_target.covs), Nt, 3, 3);
     S_t = permute(S_t, [3, 2, 1]);
     mu_t = gmm_target.means;

     % Get joint weights
     w_ts =  reshape(repmat(w_t', Ns, 1), [], 1) .* repmat(w_s, Nt, 1);

     % Compute normalization term
     w_ss = reshape(repmat(w_s', Ns, 1), [], 1) .* repmat(w_s, Ns, 1);
     Et = -corr_grad_hess([0 0 0 0 0 0], Ns, Ns, w_ss, S_s, S_s, mu_s, mu_s, 1);


     % Optimization options
     unc_options = optimoptions(@fminunc, 'Algorithm', 'trust-region', ...
                           'SpecifyObjectiveGradient',true, ...
                           'HessianFcn', 'objective',...
                           'MaxFunEvals', 1000, ...
                           'Display', 'none',...
                           'FunctionTolerance',1e-6, ...
                           'OptimalityTolerance', 1e-6,...
                           'StepTolerance', 1e-6);

     % Run optimization
     [x_opt, score, eflag, output, grad, Hess] = ...
    fminunc(@(p)corr_grad_hess(p, Nt, Ns, w_ts, S_t, S_s, mu_t, mu_s, Et),...
                            x_init, unc_options);
end

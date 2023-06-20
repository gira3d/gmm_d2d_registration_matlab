function [out] = pose_inverse(in)
  R = in(1:3, 1:3);
  T = in(1:3, 4);
  T = -1*inv(R)*T;
  R = inv(R);
  out(1:3,1:3) = R;
  out(1:3, 4) = T;
end

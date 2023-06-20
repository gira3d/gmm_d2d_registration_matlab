function [out] = pose_compose(T1, T2)
  R1 = T1(1:3, 1:3);
  t1 = T1(1:3, 4);
  R2 = T2(1:3, 1:3);
  t2 = T2(1:3, 4);

  out(1:3, 4) = t1 + R1*t2;
  out(1:3, 1:3) = R1*R2;
end

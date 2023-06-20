function [x_opt, score, grad, Hess] = isoplanar_hybrid_registration(source_file, target_file, x_init)
  [x_opt, score] = isoplanar_registration(source_file, target_file, x_init);
  [x_opt, score] = anisotropic_registration(source_file, target_file, x_opt);
end

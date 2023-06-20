function [x_opt, score, grad, Hess] = anisotropic_registration(source_file, target_file, x_init)
  gmm_target = GMM3();
  gmm_target.load(target_file);

  gmm_source = GMM3();
  gmm_source.load(source_file);

  [x_opt, score, grad, Hess] = gmm_registration(gmm_source, gmm_target, x_init);
end

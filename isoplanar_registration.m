function [x_opt, score, grad, Hess] = isoplanar_registration(source_file, target_file, x_init)
  gmm_target = GMM3();
  gmm_target.load(target_file);
  gmm_target.isoplanarCovariances()

  gmm_source = GMM3();
  gmm_source.load(source_file);
  gmm_source.isoplanarCovariances()

  [x_opt, score, grad, Hess] = gmm_registration(gmm_source, gmm_target, x_init);
end

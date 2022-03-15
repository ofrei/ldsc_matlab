function [loglike, w] = mixer_loglike(z1, z2, ell, w, N1, N2, params)
  % params.sig2_beta, params.sig2_zero, params.rho_beta, params.rho_zero

  defvec = w > 0;
  z1=z1(defvec);  z2=z2(defvec);  ell=ell(defvec);  weights=w(defvec);  N1=N1(defvec);  N2=N2(defvec);

  cov_zero = sqrt(prod(params.sig2_zero)) * params.rho_zero;
  cov_beta = sqrt(prod(params.sig2_beta)) * params.rho_beta;

  sig1 = sqrt(params.sig2_zero(1) + N1 .* ell * params.sig2_beta(1));
  sig2 = sqrt(params.sig2_zero(2) + N2 .* ell * params.sig2_beta(2));
  rho = (cov_zero + sqrt(N1 .* N2) .* ell * cov_beta) ./ (sig1 .* sig2);

  % https://mathworld.wolfram.com/BivariateNormalDistribution.html
  z = (z1.^2 ./ sig1.^2) - 2 * rho .* z1 .* z2 ./ (sig1 .* sig2) + (z2.^2 ./ sig2.^2);
  pdf = 1 ./ (2 * pi .* sig1 .* sig2 .* sqrt(1-rho.^2)) .* exp(-z ./ (2 * (1-rho.^2)));

  loglike = -sum(weights .* log(pdf));
end

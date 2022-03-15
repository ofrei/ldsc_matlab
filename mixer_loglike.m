function loglike = mixer_loglike(z1, ell, w, N1, params)
  % params.sig2_beta, params.sig2_zero

  defvec = w > 0;
  z1=z1(defvec);  ell=ell(defvec);  weights=w(defvec);  N1=N1(defvec);

  sig1 = sqrt(params.sig2_zero + N1 .* ell * params.sig2_beta);

  pdf = normpdf(z1, 0, sig1);

  loglike = -sum(weights .* log(pdf));
end

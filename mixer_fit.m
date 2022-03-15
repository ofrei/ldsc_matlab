function results = mixer_fit(z1, z2, ell, LDmat, N1, N2, M)
  % params.sig2_beta, params.sig2_zero, params.rho_beta, params.rho_zero

  w = mixer_randprune(z1, z2, ell, LDmat, N1, N2);
  defvec = w > 0;
  z1=z1(defvec);  z2=z2(defvec);  ell=ell(defvec); N1=N1(defvec);  N2=N2(defvec);
  w = w(defvec);

  sig2_beta_grid = logspace(-8, -6, 40);
  sig2_zero_grid = linspace(0.9, 1.2, 40);
  loglike1 = nan(length(sig2_beta_grid), length(sig2_zero_grid));
  loglike2 = nan(size(loglike1));

  for sig2_beta_i=1:length(sig2_beta_grid)
    fprintf('%i of %i\n', sig2_beta_i, length(sig2_beta_grid));
    for sig2_zero_i=1:length(sig2_zero_grid)
      sig2_beta = sig2_beta_grid(sig2_beta_i);
      sig2_zero = sig2_zero_grid(sig2_zero_i);
      params = struct('sig2_beta', sig2_beta, 'sig2_zero', sig2_zero);
      loglike1(sig2_beta_i, sig2_zero_i) = mixer_loglike(z1, ell, w, N1, params);
      loglike2(sig2_beta_i, sig2_zero_i) = mixer_loglike(z2, ell, w, N2, params);
    end
  end
  %imagesc(sig2_beta_grid, sig2_zero_grid, loglike1)
  %imagesc(sig2_beta_grid, sig2_zero_grid, loglike2)

  [i1, j1] = find(loglike1==min(loglike1(:)));
  [i2, j2] = find(loglike2==min(loglike2(:)));
  sig2_beta = sig2_beta_grid([i1, i2]);
  sig2_zero = sig2_zero_grid([j1, j2]);

  rho_grid = linspace(-0.99, 0.99, 40);
  loglike = nan(length(rho_grid), length(rho_grid));
  for rho_beta_i=1:length(rho_grid)
    fprintf('%i of %i\n', rho_beta_i, length(rho_grid));
    for rho_zero_i=1:length(rho_grid)
      rho_beta = rho_grid(rho_beta_i);
      rho_zero = rho_grid(rho_zero_i);
      params = struct('sig2_beta', sig2_beta, 'sig2_zero', sig2_zero, 'rho_beta', rho_beta, 'rho_zero', rho_zero);
      loglike(rho_beta_i, rho_zero_i) = mixer_loglike2(z1, z2, ell, w, N1, N2, params);
    end
  end  
  %imagesc(rho_grid, rho_grid, loglike);
  [i, j] = find(loglike==min(loglike(:)));
  rho_beta = rho_grid(i);
  rho_zero = rho_grid(j);
  
  results = struct();
  results.params = struct('sig2_beta', sig2_beta, 'sig2_zero', sig2_zero, 'rho_beta', rho_beta, 'rho_zero', rho_zero); 
  results.loglike1 = loglike1;
  results.loglike2 = loglike2;
  results.loglike = loglike;
  results.sig2_beta_grid = sig2_beta_grid;
  results.sig2_zero_grid = sig2_zero_grid;
  results.rho_grid = rho_grid;
  results.h2 = results.params.sig2_beta * double(M);
end

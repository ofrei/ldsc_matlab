function w = mixer_randprune(z1, z2, ell, LDmat, N1, N2)
  w = 0; num_randprune=10;
  for iter=1:num_randprune
    defvec = isfinite(z1+z2+ell+N1+N2);
    logpvec=rand(size(defvec));
    logpvec(~defvec) = nan;
    logpvec=FastPrune(logpvec, LDmat);
    defvec(~isfinite(logpvec)) = 0;
    w = w + defvec;
  end
  w = w / num_randprune;
end

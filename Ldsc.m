function [est, hsq] = Ldsc(self, y, x, w, N, M, step1_ii, update_func, update_weights)
    %% Matlab port of LD_Score_Regression
    % https://github.com/bulik/ldsc/blob/master/ldscore/regressions.py#L140
    %
    % Parameters
    % ----------
    % y - matrix with shape [n_snp, 1], chi square statistic for each SNP
    % x - matrix with shape [n_snp, p], total LD scores for regression
    % w - matrix with shape [n_snp, 1], LD scores computed with sum r^2
    %     taken over only those SNPs included in the regression
    % N - matrix of ints > 0 with shape [n_snp, 1],
    %     Number of individuals sampled for each SNP
    % M - float > 0, Number of SNPs used for estimating LD score (need not
    %     equal number of SNPs included in the regression)
    % update_func - lambda, update function for IRWLS
    % update_weights - lambda, update weights (used to calc. initial weight)
    %
    % self - struct with optional parameters:
    % self.intercept - float, value to constrain intercept.
    %                  NaN or missing field => unconstrained LD score regression
    % self.two_step - float, test statistic bound for use with the two-step estimator.
    %                 Not compatible with opts.intercept.
    %                 Not compatible with partitioned LD scores.
    %
    % Returns
    % -------
    % est - matrix with shape [p + 1, 1], regression weights; last element is an intercept
    % hsq - float, a heritability estimate when regression is done on a
    %       single LD score; otherwise NaN.
    %       
    % Note: this code skips the jackknife part.
    % This imply that the estimated regression coefficients are subject to
    % attenuation bias.

    if ~isfield(self, 'intercept'), self.intercept = NaN; end;
    if ~isfield(self, 'two_step'), self.two_step = NaN; end;
    if ~isfield(self, 'chisq_max'), self.chisq_max = NaN; end;
    constrain_intercept = isfinite(self.intercept);
    n_annot = size(x, 2);

    if isfinite(self.two_step) && constrain_intercept, error('twostep is not compatible with constrain_intercept'); end;
    if isfinite(self.two_step) && n_annot > 1, error('twostep not compatible with partitioned LD Score yet'); end;

    if length(N) == 1, N = repmat(N, size(y)); end;
    
    defvec = isfinite(y + sum(x, 2) + w);
    if isfinite(self.chisq_max)
        %fprintf('Removed %i SNPs with chi^2 > %.2f\n', sum(y >= chisq_max), chisq_max);
        defvec = defvec & (y < self.chisq_max);
    end
    
    y = y(defvec); x = x(defvec, :); w = w(defvec); N = N(defvec);
    if ~isempty(step1_ii), step1_ii = step1_ii(defvec); end;
    self.defvec = defvec;
    
    M_tot = sum(M);    % sum across annotation categories
    x_tot = sum(x, 2);
    tot_agg = aggregate(self, y, x_tot, N, M_tot);
    initial_w = update_weights(self, x_tot, w, N, M_tot, tot_agg, self.intercept);
    Nbar = mean(N);
    x = (repmat(N, [1 n_annot]) .* x) / Nbar;
    
    if ~constrain_intercept
        x = append_intercept(x);
        x_tot = append_intercept(x_tot);
        yp = y;
    else
        yp = y - self.intercept;
    end
    clear('y');
    
    if isfinite(self.two_step)
        yp1 = yp(step1_ii); x1 = x(step1_ii, :); w1 = w(step1_ii);
        N1 = N(step1_ii); initial_w1 = initial_w(step1_ii);
        step1_est = IRWLS(x1, yp1, @(a)update_func(self, a, x1, w1, N1, M_tot, Nbar, self.intercept, step1_ii), initial_w1);
        step1_int = step1_est(end);
        yp = yp - step1_int;
        x = remove_intercept(x);
        x_tot = remove_intercept(x_tot);
        est = IRWLS(x, yp, @(a)update_func(self, a, x_tot, w, N, M_tot, Nbar, step1_int), initial_w);
        est = [est step1_int];
        hsq = M_tot * est(1) / Nbar;
    elseif n_annot > 1  % use old_weights=True
        w = sqrt(initial_w);
        est = wls(x, yp, w);
        hsq = nan;
    else
        est = IRWLS(x, yp, @(a)update_func(self, a, x_tot, w, N, M_tot, Nbar, self.intercept), initial_w);
        hsq = M_tot * est(1) / Nbar;
    end
end

function tot_agg = aggregate(self, y, x, N, M)
    %% Calculate an initial approximation for heritability (educated guess)
    if isfinite(self.intercept), intercept = self.intercept; else intercept = self.null_intercept; end;
    num = M * (mean(y) - intercept);
    denom = mean(x .* N);
    tot_agg = num / denom;
end

function x_new = append_intercept(x)
    %% Appends an intercept term to the design matrix for a linear regression.
    n_row = size(x, 1);
    intercept = ones(n_row, 1);
    x_new = [x, intercept];
end

function x_new = remove_intercept(x)
    %% Removes the last column.
    x_new = x(:, 1:end-1);
end

function est = IRWLS(x, y, update_func, w)
    %% Iteratively re-weighted least squares (IRWLS).
    %
    % Parameters
    % ----------
    % x - matrix with shape (n, p), independent variable
    % y - matrix with shape (n, 1), dependent variable
    % update_func - transforms output of mldivide to new weights.
    % w - initial regression weights
    %
    % Returns
    % -------
    % est - matrix with shape (1, p), regression weights
    w = sqrt(w);
    for i=1:2  % update this later
        w = sqrt(update_func(wls(x, y, w)));
    end
    
    est = wls(x, y, w);
end

function est = wls(x, y, w)
    %% Weighted least squares
    w = w / sum(w);
    x = x .* repmat(w, [1, size(x, 2)]);
    y = y .* w;
    est = x \ y;  % mldivide
end
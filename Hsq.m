function [est, hsq] = Hsq(y, x, w, N, M, opts)
    %% Matlab port of Hsq estimate in LD score regression
    % https://github.com/bulik/ldsc/blob/master/ldscore/regressions.py#L336
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
    % opts - struct with optional parameters:
    % opts.intercept - float, value to constrain intercept.
    %                  NaN or missing field => unconstrained LD score regression
    % opts.two_step - float, test statistic bound for use with the two-step estimator.
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

    if ~exist('opts', 'var'), self = struct(); else self = opts; clear('opts'); end;
    
    if ~isfield(self, 'intercept'), self.intercept = NaN; end;
    if ~isfield(self, 'two_step'), self.two_step = NaN; end;
    self.constrain_intercept = isfinite(self.intercept);
    self.null_intercept = 1.0;
    self.n_annot = size(x, 2);

    % Default setting in LDSC is to use two-step estimator with cutoff 30.
    if (self.n_annot == 1) && ~self.constrain_intercept && ~isfinite(self.two_step), self.two_step = 30; end;
    
    if isfinite(self.two_step) && self.constrain_intercept, error('twostep is not compatible with constrain_intercept'); end;
    if isfinite(self.two_step) && self.n_annot > 1, error('twostep not compatible with partitioned LD Score yet'); end;
    
    if length(N) == 1, N = repmat(N, size(y)); end;
    
    defvec = isfinite(y + sum(x, 2) + w);
    if self.n_annot > 1
        chisq_max = max(0.001*max(N), 80);
        %fprintf('Removed %i SNPs with chi^2 > %.2f\n', sum(y >= chisq_max), chisq_max);
        defvec = defvec & (y < chisq_max);
    end
    
    y = y(defvec); x = x(defvec, :); w = w(defvec); N = N(defvec);
    
    if isfinite(self.two_step), step1_ii = y < self.two_step; end;
    
    M_tot = sum(M);    % sum across annotation categories
    x_tot = sum(x, 2);
    tot_agg = aggregate(self, y, x_tot, N, M_tot);
    initial_w = update_weights(self, x_tot, w, N, M_tot, tot_agg, self.intercept);
    Nbar = mean(N);
    x = (repmat(N, [1 self.n_annot]) .* x) / Nbar;
    
    if ~self.constrain_intercept
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
        step1_est = IRWLS(x1, yp1, @(a)update_func(self, a, x1, w1, N1, M_tot, Nbar), initial_w1);
        step1_int = step1_est(end);
        yp = yp - step1_int;
        x = remove_intercept(x);
        x_tot = remove_intercept(x_tot);
        est = IRWLS(x, yp, @(a)update_func(self, a, x_tot, w, N, M_tot, Nbar, step1_int), initial_w);
        hsq = M_tot * est(1) / Nbar;
    elseif self.n_annot > 1  % use old_weights=True
        w = sqrt(initial_w);
        est = wls(x, yp, w);
        hsq = nan;
    else
        est = IRWLS(x, yp, @(a)update_func(self, a, x_tot, w, N, M_tot, Nbar), initial_w);
        hsq = M_tot * est(1) / Nbar;
    end
end

function w = update_weights(self, ld, w_ld, N, M, hsq, intercept)
    %% Regression weights
    %
    % Parameters
    % ----------
    % ld - matrix with shape [n_snp, 1], LD scores
    % w_ld - matrix with shape [n_snp, 1], LD scores computed with sum r^2
    %        taken over only those SNPs included in the regression
    % N - matrix of ints > 0 with shape [n_snp, 1],
    %     Number of individuals sampled for each SNP
    % M - float > 0, Number of SNPs used for estimating LD score (need not
    %     equal number of SNPs included in the regression)
    % hsq - float in [0, 1], Heritability estimate
    %
    % Returns
    % -------
    % w - matrix with shape [n_snp, 1], Regression weights.
    %     Approx equal to reciprocal of conditional variance function.

    if ~isfinite(intercept), intercept = self.null_intercept; end;

    hsq = max(hsq, 0.0);
    hsq = min(hsq, 1.0);
    ld = max(ld, 1.0);
    w_ld = max(w_ld, 1.0);
    c = hsq * N / M;
    het_w = 1.0 ./ (2 * (intercept + c .* ld).^2);
    oc_w = 1.0 ./ w_ld;
    w = het_w .* oc_w;
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

function w = update_func(self, x, ref_ld_tot, w_ld, N, M, Nbar, intercept)
    %% Update function for IRWLS
    % x is the output of mldivide, e.g. regression coefficients
    % size(x) is (# of dimensions, 1)
    % the lst element of x is the intercept
    %
    % ~isfinite(intercept) --> free intercept
    % isfinite(intercept) --> constrained intercept
    
    if ~exist('intercept', 'var'), intercept = self.intercept; end;
    hsq = M * x(1) / Nbar;
    if ~isfinite(intercept)
        intercept = max(x(end));  % divide by zero error if intercept < 0
    else
        if size(ref_ld_tot, 2) > 1
            error('Design matrix has intercept column for constrained intercept regression!');
        end
    end

    ld = ref_ld_tot(:, 1);
    w = update_weights(self, ld, w_ld, N, M, hsq, intercept);
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
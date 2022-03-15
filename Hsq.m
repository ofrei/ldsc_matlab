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
    % opts.disable_update_weights   - true/false, whether to disable update weights
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
    %
    % TBD: chisq_max works slightly different from the original ldsc.py

    if ~exist('opts', 'var'), self = struct(); else self = opts; clear('opts'); end;
    if ~isfield(self, 'intercept'), self.intercept = NaN; end;
    if ~isfield(self, 'two_step'), self.two_step = NaN; end;
    if ~isfield(self, 'chisq_max'), self.chisq_max = NaN; end;
    if ~isfield(self, 'disable_update_weights'), self.disable_update_weights = false; end;
    self.null_intercept = 1.0;
    n_annot = size(x, 2);
    M = double(M);

    % Default setting in LDSC is to use two-step estimator with cutoff 30.
    if (n_annot > 1) && ~isfinite(self.chisq_max), self.chisq_max = max(0.001*max(N), 80); end;
    if (n_annot == 1) && ~isfinite(self.intercept) && ~isfinite(self.two_step), self.two_step = 30; end;

    if isfinite(self.two_step), step1_ii = y < self.two_step; else step1_ii = []; end;
    [est, hsq] = Ldsc(self, y, x, w, N, M, step1_ii, @update_func, @update_weights);
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

    if self.disable_update_weights
        w = w_ld;
        return
    end

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

function w = update_func(self, x, ref_ld_tot, w_ld, N, M, Nbar, intercept, ii)
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

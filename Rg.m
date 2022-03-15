function [est, rg] = Rg(z1, z2, x, w, N1, N2, M, opts)
    %% Matlab port of RG estimate in LD score regression
    % https://github.com/bulik/ldsc/blob/master/ldscore/regressions.py#L680
    %
    % Parameters
    % ----------
    % z1, z2 - matrices with shape [n_snp, 1], z statistic for each SNP (two phenotypes)
    % x - matrix with shape [n_snp, p], total LD scores for regression
    % w - matrix with shape [n_snp, 1], LD scores computed with sum r^2
    %     taken over only those SNPs included in the regression
    % N1, N2 - matrices of ints > 0 with shape [n_snp, 1],
    %     Number of individuals sampled for each SNP (two studies)
    % M - float > 0, Number of SNPs used for estimating LD score (need not
    %     equal number of SNPs included in the regression)
    % opts - struct with optional parameters:
    % opts.intercept_hsq1   - float, value to constrain intercept (first phenotype)
    % opts.intercept_hsq2   - float, value to constrain intercept (second phenotype)
    % opts.intercept_gencov - float, value to constrain intercept (covariance).
    %                         NaN or missing field => unconstrained LD score regression
    % opts.disable_update_weights   - true/false, whether to disable update weights
    % opts.two_step - float, test statistic bound for use with the two-step estimator.
    %                 Not compatible with opts.intercept_gencov.
    %                 Not compatible with partitioned LD scores.
    %
    % Returns
    % -------
    % est - matrix with shape [p + 1, 1], regression weights; last element is an intercept
    % rg  - float, a genetic corelation estimate
    %       
    % Note: this code skips the jackknife part.
    % This imply that the estimated regression coefficients are subject to
    % attenuation bias.
    %
    % TBD: chisq_max works slightly different from the original ldsc.py

    %defvec = isfinite(z1+z2+N1+N2) ; % & (x > 20);
    %z1(~defvec)= NaN;    z2(~defvec)= NaN;    N1(~defvec)= NaN;    N2(~defvec)= NaN;

    if ~exist('opts', 'var'), opts = struct(); end;
    if ~isfield(opts, 'intercept_hsq1'), opts.intercept_hsq1 = NaN; end;
    if ~isfield(opts, 'intercept_hsq2'), opts.intercept_hsq2 = NaN; end;
    if ~isfield(opts, 'intercept_gencov'), opts.intercept_gencov = NaN; end;
    if ~isfield(opts, 'two_step'), opts.two_step = NaN; end;
    if ~isfield(opts, 'chisq_max'), opts.chisq_max = NaN; end;
    if ~isfield(opts, 'disable_update_weights'), opts.disable_update_weights = false; end;
    n_annot = size(x, 2);
    M = double(M);
    
    
    opts.intercept = opts.intercept_hsq1;
    [est1, hsq1] = Hsq(z1.^2, x, w, N1, M, opts);
    opts.intercept_hsq1 = est1(end);
    opts.hsq1 = hsq1;
    fprintf('trait1: h2=%.3f intercept=%.3f\n', opts.hsq1,  opts.intercept_hsq1);
    
    opts.intercept = opts.intercept_hsq2;
    [est2, hsq2] = Hsq(z2.^2, x, w, N2, M, opts);
    opts.intercept_hsq2 = est2(end);
    opts.hsq2 = hsq2;
    fprintf('trait2: h2=%.3f intercept=%.3f\n', opts.hsq2,  opts.intercept_hsq2);
    
    opts.intercept = opts.intercept_gencov;
    [est, gencov] = Gencov(z1, z2, x, w, N1, N2, M, opts);
    rg = gencov ./ sqrt(opts.hsq1 * opts.hsq2);
    fprintf('rg=%.3f rg_intercept=%.3f\n', rg, est(end));
end

function [est, gencov] = Gencov(z1, z2, x, w, N1, N2, M, opts)
    self = opts; clear('opts');
    self.null_intercept = 1.0;
    self.N1 = N1;
    self.N2 = N2;

    n_annot = size(x, 2);
    if (n_annot == 1) && ~isfinite(self.intercept) && ~isfinite(self.two_step), self.two_step = 30; end;

    step1_ii = [];
    if isfinite(self.two_step),
        step1_ii = (z1.^2 < self.two_step & z2 .^ 2 < self.two_step);
    end;
            
    [est, gencov] = Ldsc(self, z1 .* z2, x, w, sqrt(N1 .* N2), M, step1_ii, @update_func, @update_weights);
end

function w = update_weights(self, ld, w_ld, sqrt_n1n2, M, rho_g, intercept_gencov, ii)
    %% Regression weights
    if self.disable_update_weights
        w = w_ld;
        return
    end

    if ~isfinite(intercept_gencov), intercept_gencov = self.null_intercept; end;
    if ~exist('ii', 'var'), ii=[]; end;

    h1 = max(self.hsq1, 0.0); h2 = max(self.hsq2, 0.0);
    h1 = min(     h1,   1.0); h2 = min(     h2,   1.0);
    rho_g = min(rho_g, 1.0);
    rho_g = max(rho_g, -1.0);
    ld = max(ld, 1.0);
    w_ld = max(w_ld, 1.0);
    
    N1 = self.N1; N2 = self.N2;
    if isfield(self, 'defvec'), N1 = N1(self.defvec); N2 = N2(self.defvec); end;
    if ~isempty(ii), N1 = N1(ii); N2 = N2(ii); end;
            
    a = h1 * N1 .* ld / M + self.intercept_hsq1;
    b = h2 * N2 .* ld / M + self.intercept_hsq2;
    sqrt_n1n2 = sqrt(N1 .* N2);

    c = (rho_g * sqrt_n1n2 .* ld) / M + intercept_gencov;
    het_w = 1.0 ./ (a .* b + c .^ 2);
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
    if ~exist('ii', 'var'), ii = []; end;
    rho_g = M * x(1) / Nbar;
    if ~isfinite(intercept)
        intercept = max(x(end));  % divide by zero error if intercept < 0
    end

    ld = ref_ld_tot(:, 1);
    w = update_weights(self, ld, w_ld, N, M, rho_g, intercept, ii);
end

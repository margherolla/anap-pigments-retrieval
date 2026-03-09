function out = function_anap_estimation_from_cp660(lambda, cp_660, coeffs, opts)
% FUNCTION_ANAP_ESTIMATION_FROM_CP660
%
% Generic function to derive aNAP(lambda) from cp(660), optionally
% propagating uncertainty.
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% lambda   : 1xL or Lx1 vector of wavelengths [nm]
% cp_660   : nSta x 1 vector of cp at 660 nm
% coeffs   : structure containing regression coefficients
%            Required fields:
%              coeffs.a_MC
%              coeffs.b_MC
%              coeffs.alpha1_mc
%              coeffs.beta1_mc
%              coeffs.alpha2_mc
%              coeffs.beta2_mc
%
% opts     : structure with options (all optional)
%
%   opts.use_uncertainty        = true/false
%   opts.use_cp660_uncertainty  = true/false
%   opts.cp_660_sd              = nSta x 1 vector
%   opts.lambda_ref             = 400
%   opts.lambda_s1              = 518
%   opts.lambda_s2              = 533
%   opts.random_seed            = 42
%   opts.verbose                = true
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% out : structure containing outputs
%
% Common fields:
%   out.pathway
%   out.lambda
%   out.cp_660
%   out.aNAP
%   out.aNAP400
%   out.S
%
% If uncertainty is enabled:
%   out.aNAP_median
%   out.aNAP_std
%   out.aNAP_p01
%   out.aNAP_p05
%   out.aNAP_p10
%   out.aNAP_p20
%   out.aNAP_p50
%   out.aNAP_p95
%   out.aNAP400_median
%   out.aNAP400_std
%   out.aNAP400_sem
%   out.S_median
%   out.S_std
%   out.S_sem
%   out.aNAP_all
%
% -------------------------------------------------------------------------

    %% --------------------------------------------------------------------
    % Validate inputs and defaults
    % ---------------------------------------------------------------------
    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    lambda = lambda(:)';
    cp_660 = cp_660(:);

    nSta = numel(cp_660);
    L    = numel(lambda);

    opts = set_default(opts, 'use_uncertainty', false);
    opts = set_default(opts, 'use_cp660_uncertainty', false);
    opts = set_default(opts, 'cp_660_sd', []);
    opts = set_default(opts, 'lambda_ref', 400);
    opts = set_default(opts, 'lambda_s1', 518);
    opts = set_default(opts, 'lambda_s2', 533);
    opts = set_default(opts, 'random_seed', 42);
    opts = set_default(opts, 'verbose', true);

    if opts.use_cp660_uncertainty && ~opts.use_uncertainty
        warning('opts.use_cp660_uncertainty=true but opts.use_uncertainty=false. cp(660) uncertainty will be ignored.');
    end

    if opts.use_cp660_uncertainty
        if isempty(opts.cp_660_sd)
            error('opts.cp_660_sd must be provided when opts.use_cp660_uncertainty = true.');
        end
        opts.cp_660_sd = opts.cp_660_sd(:);
        if numel(opts.cp_660_sd) ~= nSta
            error('opts.cp_660_sd must have the same length as cp_660.');
        end
    end

    required_fields = {'a_MC','b_MC','alpha1_mc','beta1_mc','alpha2_mc','beta2_mc'};
    for k = 1:numel(required_fields)
        if ~isfield(coeffs, required_fields{k})
            error('Missing coefficient field: coeffs.%s', required_fields{k});
        end
    end

    ok_cp = isfinite(cp_660) & (cp_660 > 0);
    if ~any(ok_cp)
        error('No valid positive cp(660) values found.');
    end

    dlam_ref = lambda - opts.lambda_ref;
    dLam_fit = opts.lambda_s2 - opts.lambda_s1;

    %% --------------------------------------------------------------------
    % Deterministic estimate from median coefficients
    % ---------------------------------------------------------------------
    a400_intercept = nanmedian(coeffs.a_MC(:));
    a400_slope     = nanmedian(coeffs.b_MC(:));

    a1 = nanmedian(coeffs.alpha1_mc(:));
    b1 = nanmedian(coeffs.beta1_mc(:));
    a2 = nanmedian(coeffs.alpha2_mc(:));
    b2 = nanmedian(coeffs.beta2_mc(:));

    ln_a400_det = nan(nSta,1);
    ln_a400_det(ok_cp) = a400_intercept + a400_slope .* log(cp_660(ok_cp));

    Delta_a = a1 - a2;
    Delta_b = b1 - b2;

    S_det = nan(nSta,1);
    S_det(ok_cp) = (Delta_a + Delta_b .* log(cp_660(ok_cp))) / dLam_fit;

    ln_a_det = ln_a400_det + (-S_det) .* dlam_ref;   % nSta x L
    aNAP_det = exp(ln_a_det);
    aNAP400_det = exp(ln_a400_det);

    %% --------------------------------------------------------------------
    % Initialize output with deterministic fields
    % ---------------------------------------------------------------------
    out = struct();
    out.lambda  = lambda;
    out.cp_660  = cp_660;
    out.options = opts;

    % common aliases
    out.aNAP    = aNAP_det;
    out.aNAP400 = aNAP400_det;
    out.S       = S_det;

    %% --------------------------------------------------------------------
    % CASE 1: no uncertainty
    % ---------------------------------------------------------------------
    if ~opts.use_uncertainty
        out.pathway = 'deterministic_only';

        if opts.verbose
            fprintf('Pathway: deterministic only\n');
            fprintf('Returned outputs: aNAP, aNAP400, S\n');
        end
        return
    end

    %% --------------------------------------------------------------------
    % CASE 2 and 3: uncertainty propagation
    % ---------------------------------------------------------------------
    a400_MC = coeffs.a_MC(:);
    b400_MC = coeffs.b_MC(:);
    a1_MC   = coeffs.alpha1_mc(:);
    b1_MC   = coeffs.beta1_mc(:);
    a2_MC   = coeffs.alpha2_mc(:);
    b2_MC   = coeffs.beta2_mc(:);

    Nmc = numel(a400_MC);
    if any([numel(b400_MC), numel(a1_MC), numel(b1_MC), numel(a2_MC), numel(b2_MC)] ~= Nmc)
        error('All coefficient MC vectors must have the same length.');
    end

    rng(opts.random_seed);

    aNAP_all    = nan(nSta, L, Nmc);
    S_all       = nan(nSta, Nmc);
    ln_a400_all = nan(nSta, Nmc);
    a400_all    = nan(nSta, Nmc);

    logcp_base = nan(nSta,1);
    logcp_base(ok_cp) = log(cp_660(ok_cp));

    if opts.use_cp660_uncertainty
        out.pathway = 'uncertainty_coefficients_plus_cp660';

        cp_mu = cp_660;
        cp_sd = opts.cp_660_sd;

        sigma_log = nan(nSta,1);
        good_sd = isfinite(cp_mu) & isfinite(cp_sd) & (cp_mu > 0) & (cp_sd >= 0);
        sigma_log(good_sd) = sqrt(log(1 + (cp_sd(good_sd) ./ cp_mu(good_sd)).^2));
    else
        out.pathway = 'uncertainty_coefficients_only';
    end

    for k = 1:Nmc
        a400_k = a400_MC(k);  b400_k = b400_MC(k);
        a1_k   = a1_MC(k);    b1_k   = b1_MC(k);
        a2_k   = a2_MC(k);    b2_k   = b2_MC(k);

        logcp_mc = logcp_base;

        if opts.use_cp660_uncertainty
            valid_jitter = isfinite(logcp_base) & isfinite(sigma_log);
            logcp_mc(valid_jitter) = logcp_base(valid_jitter) + ...
                                     sigma_log(valid_jitter) .* randn(sum(valid_jitter),1);
        end

        ln_a400_k = nan(nSta,1);
        S_k       = nan(nSta,1);

        valid = isfinite(logcp_mc);

        ln_a400_k(valid) = a400_k + b400_k .* logcp_mc(valid);
        S_k(valid)       = ((a1_k - a2_k) + (b1_k - b2_k) .* logcp_mc(valid)) / dLam_fit;

        ln_a_k = ln_a400_k + (-S_k) .* dlam_ref;
        aNAP_k = exp(ln_a_k);

        aNAP_all(:,:,k)   = aNAP_k;
        S_all(:,k)        = S_k;
        ln_a400_all(:,k)  = ln_a400_k;
        a400_all(:,k)     = exp(ln_a400_k);
    end

    %% --------------------------------------------------------------------
    % Monte Carlo statistics
    % ---------------------------------------------------------------------
    aNAP_median = median(aNAP_all, 3, 'omitnan');
    aNAP_std    = std(aNAP_all, 0, 3, 'omitnan');

    aNAP_prct = prctile(aNAP_all, [1 5 10 25 50 75 95], 3);
    aNAP_p01 = aNAP_prct(:,:,1);
    aNAP_p05 = aNAP_prct(:,:,2);
    aNAP_p10 = aNAP_prct(:,:,3);
    aNAP_p25 = aNAP_prct(:,:,4);
    aNAP_p75 = aNAP_prct(:,:,6);

    aNAP400_median = median(a400_all, 2, 'omitnan');
    aNAP400_std    = std(a400_all, 0, 2, 'omitnan');

    S_median = median(S_all, 2, 'omitnan');
    S_std    = std(S_all, 0, 2, 'omitnan');

    n_eff = sum(isfinite(a400_all), 2);
    aNAP400_sem = aNAP400_std ./ sqrt(max(n_eff,1));
    S_sem       = S_std ./ sqrt(max(n_eff,1));

    %% --------------------------------------------------------------------
    % Store uncertainty outputs
    % ---------------------------------------------------------------------
    out.aNAP         = aNAP_median;
    out.aNAP400      = aNAP400_median;
    out.S            = S_median;

    out.aNAP_median  = aNAP_median;
    out.aNAP_std     = aNAP_std;
    out.aNAP_p05     = aNAP_p05;
    out.aNAP_p25     = aNAP_p25;
    out.aNAP_p75     = aNAP_p75;

    out.aNAP400_median = aNAP400_median;
    out.aNAP400_std    = aNAP400_std;
    out.aNAP400_sem    = aNAP400_sem;

    out.S_median = S_median;
    out.S_std    = S_std;
    out.S_sem    = S_sem;

    out.aNAP_all = aNAP_all;
    out.aNAP400_all = a400_all;
    out.S_all = S_all;

    if opts.verbose
        fprintf('Pathway: %s\n', out.pathway);
        fprintf('Returned outputs: median aNAP + uncertainty statistics\n');
    end
end

%% =========================================================================
% Helper function
%% =========================================================================
function s = set_default(s, fieldname, value)
    if ~isfield(s, fieldname) || isempty(s.(fieldname))
        s.(fieldname) = value;
    end
end
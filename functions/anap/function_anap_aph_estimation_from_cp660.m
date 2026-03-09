function out = function_anap_aph_estimation_from_cp660(lambda, ap, cp_660, coeffs, opts)
% DERIVE_ANAP_FROM_CP660
%
% Generic function to derive aNAP(lambda) from cp(660), optionally propagate
% uncertainty, and adjust aNAP so that aph = ap - aNAP remains non-negative.
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% lambda   : 1xL or Lx1 vector of wavelengths [nm]
% ap       : nSta x L matrix of particulate absorption spectra
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
%   opts.cp_660_sd             = nSta x 1 vector (required if use_cp660_uncertainty=true)
%   opts.enforce_win1          = [410 690]
%   opts.enforce_win2          = [410 690]
%   opts.eps_pos               = 1e-4
%   opts.step_factor           = 0.95
%   opts.res_factor_min        = 0.292268426617641
%   opts.lambda_ref            = 400
%   opts.lambda_s1             = 518
%   opts.lambda_s2             = 533
%   opts.random_seed           = 42
%   opts.verbose               = true
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% out : structure containing outputs
%
% Common fields:
%   out.lambda
%   out.ap
%   out.cp_660
%   out.aNAP_base
%   out.aNAP_final
%   out.aph_final
%   out.aNAP400_base
%   out.aNAP400_final
%   out.S_used
%   out.used_step1
%   out.used_step2
%   out.needs_fallback
%   out.pathway
%
% If uncertainty is enabled:
%   out.aNAP_median
%   out.aNAP_std
%   out.aNAP_p01
%   out.aNAP_p05
%   out.aNAP_p10
%   out.aNAP_p20
%   out.aNAP400_median
%   out.aNAP400_std
%   out.aNAP400_sem
%   out.S_median
%   out.S_sem
%   out.aNAP_all           (optional large array; kept here because useful)
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
% Pathways:
%   1) No uncertainty:
%      deterministic aNAP + Step 2 fallback only
%
%   2) Uncertainty + cp(660) uncertainty:
%      MC coefficients + MC cp jitter + Step 1 + Step 2
%
%   3) Uncertainty but no cp(660) uncertainty:
%      MC coefficients only + Step 1 + Step 2
%
% -------------------------------------------------------------------------

    %% -----------------------------
    % Validate and set defaults
    % -----------------------------
    if nargin < 5
        opts = struct();
    end

    lambda = lambda(:)';
    cp_660 = cp_660(:);

    [nSta, L] = size(ap);

    if numel(lambda) ~= L
        error('Length of lambda must match number of columns of ap.');
    end
    if numel(cp_660) ~= nSta
        error('Length of cp_660 must match number of rows of ap.');
    end

    % Default options
    opts = set_default(opts, 'use_uncertainty', false);
    opts = set_default(opts, 'use_cp660_uncertainty', false);
    opts = set_default(opts, 'cp_660_sd', []);
    opts = set_default(opts, 'enforce_win1', [410 690]);
    opts = set_default(opts, 'enforce_win2', [410 690]);
    opts = set_default(opts, 'eps_pos', 1e-3);
    opts = set_default(opts, 'step_factor', 0.95);
    opts = set_default(opts, 'res_factor_min', 0.292268426617641);
    opts = set_default(opts, 'lambda_ref', 400);
    opts = set_default(opts, 'lambda_s1', 518);
    opts = set_default(opts, 'lambda_s2', 533);
    opts = set_default(opts, 'random_seed', 42);
    opts = set_default(opts, 'verbose', true);

    if opts.use_cp660_uncertainty && isempty(opts.cp_660_sd)
        error('opts.cp_660_sd must be provided when opts.use_cp660_uncertainty = true.');
    end

    if opts.use_cp660_uncertainty
        opts.cp_660_sd = opts.cp_660_sd(:);
        if numel(opts.cp_660_sd) ~= nSta
            error('opts.cp_660_sd must have same length as cp_660.');
        end
    end

    required_fields = {'a_MC','b_MC','alpha1_mc','beta1_mc','alpha2_mc','beta2_mc'};
    for k = 1:numel(required_fields)
        if ~isfield(coeffs, required_fields{k})
            error('Missing coefficient field: coeffs.%s', required_fields{k});
        end
    end

    % Masks and constants
    ok_cp = isfinite(cp_660) & (cp_660 > 0);
    if ~any(ok_cp)
        error('No valid cp(660) values found.');
    end

    enforce_win1 = lambda >= opts.enforce_win1(1) & lambda <= opts.enforce_win1(2);

    wl2_max = opts.enforce_win2(2);
    idx_wl2 = find(lambda <= wl2_max, 1, 'last');
    if isempty(idx_wl2)
        idx_wl2 = L;
    end
    enforce_win2 = lambda >= opts.enforce_win2(1) & lambda <= lambda(idx_wl2);

    dlam_ref = lambda - opts.lambda_ref;
    dLam_fit = opts.lambda_s2 - opts.lambda_s1;

    %% -----------------------------
    % Deterministic aNAP from median coefficients
    % -----------------------------
    a400_intercept = nanmedian(coeffs.a_MC, 1);
    a400_slope     = nanmedian(coeffs.b_MC, 1);

    a1 = nanmedian(coeffs.alpha1_mc, 1);
    b1 = nanmedian(coeffs.beta1_mc, 1);
    a2 = nanmedian(coeffs.alpha2_mc, 1);
    b2 = nanmedian(coeffs.beta2_mc, 1);

    ln_a400_det = nan(nSta,1);
    ln_a400_det(ok_cp) = a400_intercept + a400_slope .* log(cp_660(ok_cp));

    Delta_a = a1 - a2;
    Delta_b = b1 - b2;

    S_det = nan(nSta,1);
    S_det(ok_cp) = (Delta_a + Delta_b .* log(cp_660(ok_cp))) / dLam_fit;

    ln_a_det = ln_a400_det + (-S_det) .* dlam_ref;
    aNAP_det = exp(ln_a_det);
    aNAP400_det = exp(ln_a400_det);

    %% -----------------------------
    % Initialize outputs
    % -----------------------------
    out = struct();
    out.lambda = lambda;
    out.ap = ap;
    out.cp_660 = cp_660;
    out.options = opts;

    out.used_step1 = false(nSta,1);
    out.used_step2 = false(nSta,1);
    out.needs_fallback = false(nSta,1);

    %% =============================
    % CASE 1: NO UNCERTAINTY
    % =============================
    if ~opts.use_uncertainty
        out.pathway = 'deterministic_only';

        aNAP_base = aNAP_det;
        aNAP400_base = aNAP400_det;
        S_used = S_det;

        aph0 = ap - aNAP_base;
        needs_fallback = any(aph0(:, enforce_win2) < opts.eps_pos, 2);

        aNAP_final  = aNAP_base;
        aNAP400_adj = aNAP400_base;
        used_step2  = false(nSta,1);

        for i = find(needs_fallback).'
            ap_i   = ap(i,:);
            S_i    = S_used(i);
            a400_i = aNAP400_base(i);

            if ~(isfinite(a400_i) && isfinite(S_i) && any(isfinite(ap_i)))
                continue;
            end

            factor = 1.0;
            min_factor = opts.res_factor_min;

            base_neg  = sum((ap_i(enforce_win2) - aNAP_final(i, enforce_win2)) < opts.eps_pos);
            best_neg  = base_neg;
            best_a400 = a400_i;
            best_spec = aNAP_final(i,:);

            while factor >= min_factor
                anap_trial = (a400_i * factor) .* exp(-S_i * dlam_ref);
                aph_trial  = ap_i - anap_trial;
                neg_cnt    = sum(aph_trial(enforce_win2) < opts.eps_pos);

                if neg_cnt < best_neg
                    best_neg  = neg_cnt;
                    best_a400 = a400_i * factor;
                    best_spec = anap_trial;
                end

                if neg_cnt == 0
                    aNAP400_adj(i) = a400_i * factor;
                    aNAP_final(i,:) = anap_trial;
                    used_step2(i) = true;
                    break;
                end

                factor = factor * opts.step_factor;
            end

            if ~used_step2(i) && best_neg < base_neg
                aNAP400_adj(i) = best_a400;
                aNAP_final(i,:) = best_spec;
                used_step2(i) = true;
            end
        end

        out.aNAP_base     = aNAP_base;
        out.aNAP_final    = aNAP_final;
        out.aph_final     = ap - aNAP_final;
        out.aNAP400_base  = aNAP400_base;
        out.aNAP400_final = aNAP400_adj;
        out.S_used        = S_used;
        out.needs_fallback = needs_fallback;
        out.used_step2    = used_step2;

        if opts.verbose
            fprintf('No-uncertainty path: Step 2 applied to %d/%d stations\n', nnz(used_step2), nSta);
        end

        return
    end

    %% =============================
    % CASE 2 and 3: UNCERTAINTY PATHS
    % =============================
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
            logcp_mc(valid_jitter) = logcp_base(valid_jitter) + sigma_log(valid_jitter) .* randn(sum(valid_jitter),1);
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

    % MC statistics
    aNAP_median = median(aNAP_all, 3, 'omitnan');
    aNAP_std    = std(aNAP_all, 0, 3, 'omitnan');

    aNAP_prct = prctile(aNAP_all, [1 5 25 75], 3);
    aNAP_p01 = aNAP_prct(:,:,1);
    aNAP_p05 = aNAP_prct(:,:,2);
    aNAP_p25 = aNAP_prct(:,:,3);
    aNAP_p75 = aNAP_prct(:,:,4);

    aNAP400_median = median(a400_all, 2, 'omitnan');
    aNAP400_std    = std(a400_all, 0, 2, 'omitnan');
    S_median       = median(S_all,   2, 'omitnan');

    n_eff = sum(isfinite(a400_all), 2);
    aNAP400_sem = std(a400_all, 0, 2, 'omitnan') ./ sqrt(max(n_eff,1));
    S_sem       = std(S_all,    0, 2, 'omitnan') ./ sqrt(max(n_eff,1));

    anap_before = aNAP_median;

    %% -----------------------------
    % STEP 1: p05 replacement if it solves the negativity
    % -----------------------------
    anap_step1 = aNAP_median;
    used_step1 = false(nSta,1);

    for i = 1:nSta
        api  = ap(i,:);
        ai50 = aNAP_median(i,:);
        ai05 = aNAP_p05(i,:);

        if ~any(isfinite(api(enforce_win1))) || ~any(isfinite(ai50(enforce_win1)))
            continue;
        end

        if any(ai50(enforce_win1) > api(enforce_win1))
            if all(ai05(enforce_win1) <= api(enforce_win1))
                anap_step1(i,:) = max(ai05, 0);
                used_step1(i) = true;
            end
        end
    end

    aph_after_step1 = ap - anap_step1;
    needs_fallback  = any(aph_after_step1(:, enforce_win2) < opts.eps_pos, 2);

    %% -----------------------------
    % STEP 2: multiplicative fallback on aNAP(400)
    % -----------------------------
    anap_final = anap_step1;
    anap400_adj = get_a400_from_spectra(lambda, anap_final, opts.lambda_ref);
    used_step2 = false(nSta,1);

    for i = find(needs_fallback).'
        ap_i   = ap(i,:);
        S_i    = S_median(i);
        a400_i = anap400_adj(i);
        sig_i  = aNAP400_std(i);

        if ~(isfinite(a400_i) && isfinite(S_i) && isfinite(sig_i) && any(isfinite(ap_i)))
            continue;
        end

        factor_uncert = (a400_i - sig_i) / a400_i;
        factor_uncert = max(0, min(1, factor_uncert));

        % Here I keep using the global residual minimum as actual lower bound
        min_factor = opts.res_factor_min;

        factor = 1.0;
        base_neg  = sum((ap_i(enforce_win2) - anap_final(i, enforce_win2)) < opts.eps_pos);
        best_neg  = base_neg;
        best_a400 = a400_i;
        best_spec = anap_final(i,:);

        while factor >= min_factor
            anap_trial = (a400_i * factor) .* exp(-S_i * dlam_ref);
            aph_trial  = ap_i - anap_trial;

            neg_cnt = sum(aph_trial(enforce_win2) < opts.eps_pos);

            if neg_cnt < best_neg
                best_neg  = neg_cnt;
                best_a400 = a400_i * factor;
                best_spec = anap_trial;
            end

            if neg_cnt == 0
                anap400_adj(i) = a400_i * factor;
                anap_final(i,:) = anap_trial;
                used_step2(i) = true;
                break;
            end

            factor = factor * opts.step_factor;
        end

        if ~used_step2(i) && best_neg < base_neg
            anap400_adj(i) = best_a400;
            anap_final(i,:) = best_spec;
            used_step2(i) = true;
        end
    end

    %% -----------------------------
    % Store outputs
    % -----------------------------
    out.aNAP_base      = anap_before;
    out.aNAP_final     = anap_final;
    out.aph_final      = ap - anap_final;
    out.aNAP400_base   = aNAP400_median;
    out.aNAP400_final  = anap400_adj;
    out.S_used         = S_median;

    out.aNAP_median    = aNAP_median;
    out.aNAP_std       = aNAP_std;
    out.aNAP_p01       = aNAP_p01;
    out.aNAP_p05       = aNAP_p05;
    out.aNAP_p25       = aNAP_p25;
    out.aNAP_p75       = aNAP_p75;

    out.aNAP400_median = aNAP400_median;
    out.aNAP400_std    = aNAP400_std;
    out.aNAP400_sem    = aNAP400_sem;

    out.S_median       = S_median;
    out.S_sem          = S_sem;

    out.aNAP_all       = aNAP_all;   % can be large, but useful
    out.used_step1     = used_step1;
    out.used_step2     = used_step2;
    out.needs_fallback = needs_fallback;

    if opts.verbose
        fprintf('%s\n', ['Pathway: ' out.pathway]);
        fprintf('Step 1 applied to: %d/%d stations\n', nnz(used_step1), nSta);
        fprintf('Step 2 applied to: %d/%d stations\n', nnz(used_step2), nSta);
    end
end

%% =========================================================================
% Helper functions
%% =========================================================================
function s = set_default(s, fieldname, value)
    if ~isfield(s, fieldname) || isempty(s.(fieldname))
        s.(fieldname) = value;
    end
end

function a400 = get_a400_from_spectra(lambda, spectra, lambda_ref)
    lambda = lambda(:)';
    [nSta, ~] = size(spectra);
    a400 = nan(nSta,1);

    if any(abs(lambda - lambda_ref) < 1e-10)
        idx = find(abs(lambda - lambda_ref) < 1e-10, 1, 'first');
        a400 = spectra(:, idx);
    else
        valid_rows = ~all(isnan(spectra), 2);
        tmp = interp1(lambda(:), spectra(valid_rows,:)', lambda_ref, 'linear', 'extrap')';
        a400(valid_rows) = tmp;
    end
end
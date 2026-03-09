%% RUN_EXAMPLE_FULL_WORKFLOW
% Full workflow:
% cp(660) -> aNAP -> aPH -> Gaussian decomposition -> pigment retrieval

clear; clc;

%% ------------------------------------------------------------------------
% Locate repository root
% -------------------------------------------------------------------------

repo_dir = 'C:\Users\margh\Tara_Europa\aNAP_paper_C13\5_new_algorithm\anap-pigments-retrieval';

% Add functions folder to MATLAB path
addpath(genpath(fullfile(repo_dir,'functions')))

%% ------------------------------------------------------------------------
% Load example dataset and coefficients
% -------------------------------------------------------------------------
load(fullfile(repo_dir,'example_data','example_dataset_full.mat'));
load(fullfile(repo_dir,'coefficients','coefficients_anap.mat'));

data = example_dataset_full;
coeffs = coefficients_anap;

%% ------------------------------------------------------------------------
% Build inputs
% -------------------------------------------------------------------------
lambda = data.lambda(:)';
ap     = data.ap;
cp_660 = data.cp_660;

%% ------------------------------------------------------------------------
% Recommended example case:
% uncertainty with cp(660) uncertainty
% -------------------------------------------------------------------------
opts = struct();
opts.enforce_win1 = [410 690];
opts.enforce_win2 = [410 690];
opts.eps_pos = 1e-4;
opts.step_factor = 0.95;
opts.res_factor_min = 0.292268426617641;
opts.lambda_ref = 400;
opts.lambda_s1 = 518;
opts.lambda_s2 = 533;
opts.random_seed = 42;
opts.verbose = true;

opts.use_uncertainty = true;
opts.use_cp660_uncertainty = true;

if isfield(data, 'cp_660_sd')
    opts.cp_660_sd = data.cp_660_sd;
else
    error('Field "cp_660_sd" not found in example_dataset_full.mat.');
end

%% ------------------------------------------------------------------------
% Run aNAP / aPH estimation
% -------------------------------------------------------------------------

out_anap_aph = function_anap_aph_estimation_from_cp660(lambda, ap, cp_660, coeffs, opts);

disp(' ')
disp('Full workflow - aNAP/aPH estimation finished')
disp(['Step 1 used = ' num2str(nnz(out_anap_aph.used_step1))])
disp(['Step 2 used = ' num2str(nnz(out_anap_aph.used_step2))])

%% ------------------------------------------------------------------------
% Optional diagnostic plot
% -------------------------------------------------------------------------
station_to_plot = min(179, size(ap,1));

figure('Color','w'); hold on; grid on;

plot(lambda, ap(station_to_plot,:), 'k-', 'LineWidth', 1.5, ...
    'DisplayName','a_p');

if opts.use_uncertainty
    plot(lambda, out_anap_aph.aNAP_median(station_to_plot,:), '--', ...
        'LineWidth', 1.5, 'DisplayName','a_{NAP} median');

    plot(lambda, out_anap_aph.aNAP_p05(station_to_plot,:), ':', ...
        'LineWidth', 1.5, 'DisplayName','a_{NAP} p05');
else
    plot(lambda, out_anap_aph.aNAP_base(station_to_plot,:), '--', ...
        'LineWidth', 1.5, 'DisplayName','a_{NAP} base');
end

plot(lambda, out_anap_aph.aNAP_final(station_to_plot,:), '-', ...
    'LineWidth', 1.8, 'DisplayName','a_{NAP} final');

plot(lambda, out_anap_aph.aph_final(station_to_plot,:), '-', ...
    'LineWidth', 1.5, 'DisplayName','a_{ph} final');

xlabel('\lambda (nm)');
ylabel('a(\lambda) [m^{-1}]');
title(sprintf('Station %d - %s', station_to_plot, out_anap_aph.pathway));
legend('Location','best');

%% ------------------------------------------------------------------------
% Store outputs in a total structure
% -------------------------------------------------------------------------
data_tot = struct();
data_tot.initial_dataset = data;
data_tot.anap_aph_output = out_anap_aph;

%% ========================================================================
% GAUSSIAN DECOMPOSITION FROM aPH(lambda)
% ========================================================================

aph   = out_anap_aph.aph_final;
wl    = out_anap_aph.lambda(:)';
wlunc = wl;
apunc = data.ap_sd;
acs   = 0;

[amps_aph, compspec_aph, sumspec_aph] = spectral_decomp_aph(wl, aph, wlunc, apunc, acs);

% Gaussian peak wavelengths in the same order returned by spectral_decomp_aph
peaks = [434 453 470 492 523 550 584 617 638 660 675];

if size(amps_aph,2) ~= numel(peaks)
    error('Number of columns in amps_aph (%d) does not match number of peaks (%d).', ...
        size(amps_aph,2), numel(peaks));
end

for i = 1:numel(peaks)
    data_tot.gaus_output.(sprintf('agaus_%d', peaks(i))) = amps_aph(:, i);
end

data_tot.gaus_output.amps_gaus     = amps_aph;
data_tot.gaus_output.compspec_gaus = compspec_aph;
data_tot.gaus_output.sumspec_gaus  = sumspec_aph;

disp('Gaussian decomposition finished')

%% ========================================================================
% DERIVE PIGMENT CONCENTRATIONS
% ========================================================================

table_pigments_coeff = readtable(fullfile(repo_dir, 'coefficients', 'table_pigments_coeff.csv'));

requiredCoeffCols = {'A','B'};
for i = 1:numel(requiredCoeffCols)
    if ~ismember(requiredCoeffCols{i}, table_pigments_coeff.Properties.VariableNames)
        error('Missing required column "%s" in table_pigments_coeff.csv.', requiredCoeffCols{i});
    end
end

A = table_pigments_coeff.A;
B = table_pigments_coeff.B;

stations = (1:size(aph,1))';

% 550 nm excluded because not associated with pigments used here
peaks_used = peaks(peaks ~= 550);

nSta   = size(aph,1);
nPeaks = numel(peaks_used);
agaus  = nan(nSta, nPeaks);

for i = 1:nPeaks
    field_name = sprintf('agaus_%d', peaks_used(i));

    if ~isfield(data_tot.gaus_output, field_name)
        error('Field "%s" not found in data_tot.gaus_output.', field_name);
    end

    agaus(:, i) = data_tot.gaus_output.(field_name);
end

if size(agaus,1) ~= numel(stations)
    error('Number of rows in agaus does not match number of stations.');
end

if numel(A) ~= numel(B)
    error('Coefficient vectors A and B must have the same length.');
end

pigm_der = derive_pigm_gaus(agaus, stations, A, B);

data_tot.gaus_output.pigm_der = pigm_der;

disp('Pigment derivation finished')
disp('Full workflow completed successfully')

% Optional save
% save(fullfile('example_data','example_output_full_workflow.mat'), 'data_tot');
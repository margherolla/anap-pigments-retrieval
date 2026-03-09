%% RUN_EXAMPLE_ANAP_ONLY
% Derive aNAP(lambda) from cp(660) using one selected mode:
%   1) deterministic
%   2) coefficient uncertainty only
%   3) coefficient + cp(660) uncertainty

clear; clc;

%% ------------------------------------------------------------------------
% USER CHOICE
% -------------------------------------------------------------------------
workflow_mode = 'full_uncertainty';
% Options:
%   'deterministic'
%   'coeff_uncertainty'
%   'full_uncertainty'

station_id = 180;   % station to plot (automatically limited if too large)

%% ------------------------------------------------------------------------
% Locate repository root and add functions to path
% -------------------------------------------------------------------------
repo_dir = 'C:\Users\margh\Tara_Europa\aNAP_paper_C13\5_new_algorithm\anap-pigments-retrieval\main';

repo_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_dir, 'functions')))

%% ------------------------------------------------------------------------
% Load dataset and coefficients
% -------------------------------------------------------------------------
load(fullfile(repo_dir, 'example_data', 'example_dataset_anap_only.mat'));
load(fullfile(repo_dir, 'coefficients', 'coefficients_anap.mat'));

data   = example_dataset_anap_only;
coeffs = coefficients_anap;

%% ------------------------------------------------------------------------
% Inputs
% -------------------------------------------------------------------------
lambda = data.lambda(:)';
cp_660 = data.cp_660(:);

nSta = numel(cp_660);
station_id = min(station_id, nSta);

%% ------------------------------------------------------------------------
% Common options
% -------------------------------------------------------------------------
opts = struct();
opts.lambda_ref  = 400;
opts.lambda_s1   = 518;
opts.lambda_s2   = 533;
opts.random_seed = 42;
opts.verbose     = true;

%% ------------------------------------------------------------------------
% Select workflow mode
% -------------------------------------------------------------------------
switch lower(workflow_mode)

    case 'deterministic'
        opts.use_uncertainty = false;
        opts.use_cp660_uncertainty = false;

    case 'coeff_uncertainty'
        opts.use_uncertainty = true;
        opts.use_cp660_uncertainty = false;

    case 'full_uncertainty'
        opts.use_uncertainty = true;
        opts.use_cp660_uncertainty = true;

        if isfield(data, 'cp_660_sd')
            opts.cp_660_sd = data.cp_660_sd(:);
        else
            error(['workflow_mode = "full_uncertainty" requires ', ...
                   'data.cp_660_sd in example_dataset_anap_only.mat.']);
        end

    otherwise
        error('Unknown workflow_mode: %s', workflow_mode);
end

%% ------------------------------------------------------------------------
% Run aNAP estimation
% -------------------------------------------------------------------------
out_anap = function_anap_estimation_from_cp660(lambda, cp_660, coeffs, opts);

disp(' ')
disp('aNAP estimation finished')
disp(['Selected mode: ' workflow_mode])
disp(['Returned pathway: ' out_anap.pathway])
disp(['Size of aNAP output: ' num2str(size(out_anap.aNAP,1)) ' x ' num2str(size(out_anap.aNAP,2))])

%% ------------------------------------------------------------------------
% Adaptive plot for one station
% -------------------------------------------------------------------------
figure('Color','w'); hold on; grid on;

% Main aNAP curve (always available)
plot(lambda, out_anap.aNAP(station_id,:), 'k-', 'LineWidth', 1.8, ...
    'DisplayName', 'a_{NAP}');

% Add uncertainty envelope only if available
if isfield(out_anap, 'aNAP_p05')
    plot(lambda, out_anap.aNAP_p05(station_id,:), '--', 'LineWidth', 1.5, ...
        'DisplayName', 'a_{NAP} p05');
end

if isfield(out_anap, 'aNAP_p25')
    plot(lambda, out_anap.aNAP_p25(station_id,:), ':', 'LineWidth', 1.5, ...
        'DisplayName', 'a_{NAP} p25');
end

if isfield(out_anap, 'aNAP_p75')
    plot(lambda, out_anap.aNAP_p75(station_id,:), ':', 'LineWidth', 1.5, ...
        'DisplayName', 'a_{NAP} p75');
end


xlabel('\lambda (nm)');
ylabel('a_{NAP}(\lambda) [m^{-1}]');
title(sprintf('Station %d - %s', station_id, strrep(out_anap.pathway, '_', ' ')));
legend('Location','best');


% =========================================================================
% Title: Spectral Feature Extraction from Simulated LFP/EEG Signals
%
% Author: Gabriele Mancini
% Date: 2025-12-10
%
% Description:
% This script processes simulated neural data to extract periodic and
% aperiodic spectral features from LFP or EEG proxy signals. LFP/EEG signals
% are reconstructed from synaptic currents and analyzed to quantify
% spectral properties (PSD slope, gamma power, Hurst exponent) alongside
% firing-rate statistics and network parameters.
%
% Associated manuscript:
% Bertelsen N. et al.,
% "Electrophysiologically-defined excitation–inhibition autism neurosubtypes",
% Nature Neuroscience (under review).
%
% License: MIT License
% Contact: gabriele.mancini@oiit.it
% =========================================================================

clear; close all; clc;

%% ========================================================================
% Paths and configuration
% ========================================================================

trials = 2;
root_path  = '/Users/gabrielemancini/Desktop/Bertelsen Github/In-silico';
path_sim   = fullfile(root_path, 'Simulations/Simulation results/g_and_nu/data');
path_stats = fullfile(root_path, 'Simulations/Simulation results/g_and_nu/metadata/stats_db.csv');
path_params= fullfile(root_path, 'Simulations/Simulation results/g_and_nu/metadata/parameter_db.csv');

addpath('../../functions');

% If true, analyze true LFP; otherwise compute EEG proxy from synaptic currents
LFP_true = true;

%% ========================================================================
% Load simulated LFP data
% ========================================================================

files = dir(fullfile(path_sim, '*.npy'));
files = {files.name};

LFP = cell(1, numel(files));
ID  = cell(1, numel(files));

for k = 1:numel(files)
    fprintf('\rLoading LFP file %d / %d (%.0f%%)', k, numel(files), 100*k/numel(files));
    LFP{k} = readNPY(fullfile(path_sim, files{k}));
    ID{k}  = extractAfter(erase(files{k}, '.npy'), 'ID');
end
fprintf('\n');

%% ========================================================================
% Load metadata
% ========================================================================

statsdb     = readmatrix(path_stats);     statsdb     = statsdb(2:end, :);
parameterdb = readmatrix(path_params);   parameterdb = parameterdb(2:end, 3:end);

%% ========================================================================
% Analysis parameters
% ========================================================================

fs              = 500;           % sampling frequency (Hz)
f_range         = [1, 80];
f_range_gamma   = [30, 50];
f_range_slope   = [20, 90];

params.tapers   = [3 5];         % multitaper parameters
Hparams         = get_hurst_parameters;

%% ========================================================================
% EEG proxy parameters (non-causal ERWS2 model)
% ========================================================================

ERWS2_params_noncausal = [-0.6, 0.1, -0.4, -1.9, 0.6, 3.0, 1.4, 1.7, 0.2];

%% ========================================================================
% Feature extraction loop
% ========================================================================

count = 1;

for i = 1:numel(LFP)

    id     = str2double(ID{i});
    nu_val = parameterdb(parameterdb(:,59) == id, 53);

    % Exclude high-nu regimes
    if nu_val >= 10
        continue
    end

    for trial = 1:trials

        % ---------------------------------------------------------------
        % Signal reconstruction
        % ---------------------------------------------------------------
        if LFP_true
            LFP_trace = squeeze(LFP{i}(1, trial, :));
        else
            AMPA_trace = squeeze(LFP{i}(2, trial, :));
            GABA_trace = squeeze(LFP{i}(3, trial, :));
            LFP_trace  = erws2(AMPA_trace(:), GABA_trace(:), ...
                ERWS2_params_noncausal, ...
                parameterdb(parameterdb(:,59)==id, 53));
        end

        % ---------------------------------------------------------------
        % Spectral features
        % ---------------------------------------------------------------
        [tot_pow, ~, ~] = compute_power(LFP_trace, fs, f_range);
        tot_power(trial,count) = tot_pow;

        [slope,~,~,~] = compute_slope_fooof_simulations(LFP_trace, fs, f_range_slope);
        slope_LFP_fooof(trial,count) = slope;

        gamma_pow_LFP(trial,count) = compute_psd_area(LFP_trace, fs, f_range_gamma);

        [S, ~] = mtspectrumc(LFP_trace, params);
        PSDs(trial,count,:) = S';

        % ---------------------------------------------------------------
        % Hurst exponent
        % ---------------------------------------------------------------
        H_LFP(trial,count) = bfn_mfin_ml(detrend(LFP_trace), ...
            'filter', Hparams.filter, ...
            'lb', Hparams.lb, ...
            'ub', Hparams.ub);

        % ---------------------------------------------------------------
        % Network statistics
        % ---------------------------------------------------------------
        exc  = statsdb(statsdb(:,10)==id,1);
        inh  = statsdb(statsdb(:,10)==id,5);
        isi_ = statsdb(statsdb(:,10)==id,2);
        corr_= statsdb(statsdb(:,10)==id,4);

        fr_exc(trial,count) = exc(trial);
        fr_inh(trial,count) = inh(trial);
        isi(trial,count)    = isi_(trial);
        corr_coef(trial,count) = corr_(trial);

        % ---------------------------------------------------------------
        % Model parameters
        % ---------------------------------------------------------------
        g_param(trial,count) = -parameterdb(parameterdb(:,59)==id,10) ./ ...
                                parameterdb(parameterdb(:,59)==id,13);
        nu(trial,count) = nu_val;

        count = count + 1;
    end
end

%% ========================================================================
% Average across trials
% ========================================================================

avg = @(x) (x(1,1:2:end-1) + x(2,2:2:end)) / 2;

FR_exc      = avg(fr_exc);
FR_inh      = avg(fr_inh);
isi_base    = avg(isi);
corr_base   = avg(corr_coef);
g           = avg(g_param);
nu_0        = avg(nu);
H_base      = avg(H_LFP);
gamma_base  = avg(gamma_pow_LFP);
slope_base  = avg(slope_LFP_fooof);
Pow_base    = avg(tot_power);

PSDs_mean = (PSDs(1,1:2:end-1,:) + PSDs(2,2:2:end,:)) / 2;
PSDS = squeeze(PSDs_mean)';


%% ========================================================================
% Save results
% ========================================================================

output_folder = fullfile(root_path, 'Data_to_plot');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

if LFP_true
    save(fullfile(output_folder, 'analysis_output.mat'));
else
    save(fullfile(output_folder, 'analysis_output_EEG.mat'));
end

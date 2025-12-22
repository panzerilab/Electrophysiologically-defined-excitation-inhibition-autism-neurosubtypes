function tab2output = compute_features_remove_channels(datafile, outdir, f_range, f_range_fooof)

% compute_hurst - compute the Hurst exponent
% 
% INPUT
%   datafile = full filename to a .mat file.
%   outdir = directory to save output csv files
% 
% OUTPUT
%   tab2output = a table of H values and other parameters
% 
% Example usage:
% 
%   datafile = '/Users/mlombardo/Dropbox/data/gozzi_dreadd_lfp/CAMK/ds220222a_LFP_data_chunked.mat';
%   outdir = '/Users/mlombardo/Dropbox/Manuscripts/CMI_EEG_H/gozzi_dreadd_lfp/CAMK/results';
%   tab2output = compute_hurst(datafile, outdir);
% 
% -- written by mvlombardo -- 06.12.2023

%% load data

load(datafile);
data = LFP_data_chunked;

%% PCA reduce channels


channels_idx = [1,4,8,10,12,14];

% Baseline
baseline_windows = size(data.baseline, 3);
baseline_data = nan(size(data.baseline, 1), baseline_windows);
for iwin = 1:baseline_windows
    [coeff, score, latent, tsquared, explained] = pca(squeeze(data.baseline(:,channels_idx,iwin)));
%     disp(explained(1));
    baseline_data(:, iwin) = score(:, 1); 
end % for iwin

% Transition
transition_windows = 1:20;
transition_data = nan(size(data.cno, 1), length(transition_windows));
for iwin = 1:length(transition_windows)
    window2use = transition_windows(iwin);
    [coeff, score, latent, tsquared, explained] = pca(squeeze(data.cno(:,channels_idx,window2use)));
%     disp(explained(1));
    transition_data(:, iwin) = score(:, 1); 
end % for iwin

% Treatment
treatment_windows = (transition_windows(end)+1):size(data.cno, 3);
treatment_data = nan(size(data.cno, 1), length(treatment_windows));
for iwin = 1:length(treatment_windows)
    window2use = treatment_windows(iwin);
    [coeff, score, latent, tsquared, explained] = pca(squeeze(data.cno(:,channels_idx,window2use)));
%     disp(explained(1));
    treatment_data(:, iwin) = score(:, 1); 
end % for iwin  % <-- Correctly closed the loop


%% Compute Hurst and other features

Hparams = get_hurst_parameters;

%% compute H

Hbaseline = bfn_mfin_ml(baseline_data, ...
    'filter', Hparams.filter, ...
    'lb',Hparams.lb, ...
    'ub',Hparams.ub, ...
    'verbose',Hparams.verbose);

Htransition = bfn_mfin_ml(transition_data, ...
    'filter', Hparams.filter, ...
    'lb',Hparams.lb, ...
    'ub',Hparams.ub, ...
    'verbose',Hparams.verbose);

Htreatment = bfn_mfin_ml(treatment_data, ...
    'filter', Hparams.filter, ...
    'lb',Hparams.lb, ...
    'ub',Hparams.ub, ...
    'verbose',Hparams.verbose);

% Sampling rate
fs = data.fs;
 
% Compute slope
slope_baseline = -compute_spectral_slope(baseline_data, fs, f_range);
slope_transition = -compute_spectral_slope(transition_data, fs, f_range);
slope_treatment = -compute_spectral_slope(treatment_data, fs, f_range);

% Compute FOOOF parameters for each condition
[slope_fooof_baseline, freq_baseline, amplitude_baseline, bw_baseline, n_peaks_b] = compute_spectral_slope_fooof(baseline_data, fs, f_range_fooof);
[slope_fooof_transition, freq_transition, amplitude_transition, bw_transition, n_peaks_tra] = compute_spectral_slope_fooof(transition_data, fs, f_range_fooof);
[slope_fooof_treatment, freq_treatment, amplitude_treatment, bw_treatment, n_peaks_tre] = compute_spectral_slope_fooof(treatment_data, fs, f_range_fooof);

[tot_pow_b, gamma_pow_b, low_pow_b, fx_b, PSD_b] = compute_power(baseline_data, fs, f_range);
[tot_pow_tra, gamma_pow_tra, low_pow_tra, fx_tra, PSD_tra] = compute_power(transition_data, fs, f_range);
[tot_pow_tre, gamma_pow_tre, low_pow_tre, fx_tre, PSD_tre] = compute_power(treatment_data, fs, f_range);


[delta_b, theta_b, alpha_b, beta_b, gamma_b,broad1_b,broad2_b] = compute_bands(baseline_data, fs);
[delta_tra, theta_tra, alpha_tra, beta_tra, gamma_tra,broad1_tra,broad2_tra] = compute_bands(transition_data, fs);
[delta_tre, theta_tre, alpha_tre, beta_tre, gamma_tre,broad1_tre,broad2_tre] = compute_bands(treatment_data, fs);



%% Combine and organize data

slope_fit_all = [slope_baseline'; slope_transition'; slope_treatment'];

% Consolidate all features in arrays
slope_all = [slope_fooof_baseline'; slope_fooof_transition'; slope_fooof_treatment'];
freq_all = [freq_baseline'; freq_transition'; freq_treatment'];
amplitude_all = [amplitude_baseline'; amplitude_transition'; amplitude_treatment'];
bw_all = [bw_baseline'; bw_transition'; bw_treatment'];
n_peaks_all = [n_peaks_b'; n_peaks_tra'; n_peaks_tre'];

tot_pow_all = [tot_pow_b'; tot_pow_tra'; tot_pow_tre'];
gamma_pow_all = [gamma_pow_b'; gamma_pow_tra'; gamma_pow_tre'];
low_pow_all = [low_pow_b'; low_pow_tra'; low_pow_tre'];

delta_all = [delta_b'; delta_tra'; delta_tre'];
theta_all = [theta_b'; theta_tra'; theta_tre'];
alpha_all = [alpha_b'; alpha_tra'; alpha_tre'];
beta_all = [beta_b'; beta_tra'; beta_tre'];
gamma_all = [gamma_b'; gamma_tra'; gamma_tre'];

broad1_all = [broad1_b'; broad1_tra'; broad1_tre'];
broad2_all = [broad2_b'; broad2_tra'; broad2_tre'];

fx_all_app = [fx_b; fx_tra; fx_tre];
PSD_all_app = [PSD_b; PSD_tra; PSD_tre];

fx_all = mat2cell(fx_all_app, ones(size(fx_all_app,1),1), size(fx_all_app,2));
PSD_all = mat2cell(PSD_all_app, ones(size(PSD_all_app,1),1), size(PSD_all_app,2));

Houtput = num2cell([Hbaseline';Htransition';Htreatment']); 

% Condition labels
cond_labels = [repmat({'baseline'}, length(slope_baseline), 1);
               repmat({'transition'}, length(slope_transition), 1);
               repmat({'treatment'}, length(slope_treatment), 1)];

% Generate ID, group, window, and total labels
[f_path, f_name, f_ext] = fileparts(datafile);
end_idx = strfind(f_name, '_LFP_data_chuncked');
id2use = f_name(1:(end_idx-1));
id_labels = repmat({id2use}, length(cond_labels), 1);
grp_labels = repmat({data.type}, length(cond_labels), 1);

% Window and total window labels
window_labels = num2cell([[1:length(slope_baseline)]'; 
                          [1:length(slope_transition)]'; 
                          [1:length(slope_treatment)]']);
total_windows_labels = num2cell(1:length(cond_labels))';


%% Create final table

tab2output = cell2table([num2cell(Houtput), num2cell(slope_all), num2cell(freq_all), ...
                         num2cell(amplitude_all), num2cell(bw_all), ...
                         cond_labels, id_labels, grp_labels, ...
                         window_labels, total_windows_labels, ...
                         num2cell(tot_pow_all), num2cell(gamma_pow_all), num2cell(low_pow_all), num2cell(slope_fit_all), ...
                         num2cell(delta_all), num2cell(theta_all), num2cell(alpha_all), num2cell(beta_all), num2cell(gamma_all), num2cell(broad1_all), num2cell(broad2_all), num2cell(fx_all), num2cell(PSD_all), num2cell(n_peaks_all)], ...
                        'VariableNames', {'H', 'Slope', 'Frequency', 'Amplitude', 'Bandwidth', ...
                                          'Condition', 'ID', 'Group', ...
                                          'Window_Num_Per_Cond', 'Window_Total', 'tot_pow', 'gamma_pow', 'low_pow', 'slope_fit', ...
                                          'delta', 'theta', 'alpha', 'beta', 'gamma','broad1','broad2','fx','PSD','n_peaks'});

%% Export table

% Ensure output directory exists
output_dir = fullfile(outdir, 'features', 'camk');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Define output file path
output_file = fullfile(outdir, sprintf('%s_features.csv', id2use));

% Write the table to file
writetable(tab2output, output_file);


end 

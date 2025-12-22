% =========================================================================
% Title: DATA INTEGRATION & NORMALIZATION SCRIPT
%
% Author: Gabriele Mancini
% Date: 2025-12-10
%
% Description:
%  Loads Multi-Unit Activity (MUA) and Spectral Feature (LFP) data.
%  Performs data cleaning, channel selection, Z-score normalization against
%  baseline periods, and aggregates data into subject-specific structures.
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
%  1. CONFIGURATION
%  ========================================================================

% --- Paths ---
PATH_MUA_ROOT  = '../DATA/MUA';
PATH_FEAT_ROOT = '../results/features';
PATH_OUT       = '../results/tables';

% Ensure output directory exists
if ~exist(PATH_OUT, 'dir')
    mkdir(PATH_OUT);
    fprintf('Created output directory: %s\n', PATH_OUT);
end

% --- Experiment Definitions ---
% Define groups, folders, and specific exclusions
% ID Codes: 0=HM4Di, 1=CaMKII, 2=PV

% 1. HM4Di
experiments(1).name     = 'hm4di'; 
experiments(1).folder   = 'hm4d'; % Folder name on disk
experiments(1).exclude  = {'fr190610b'}; 

% 2. CaMKII
experiments(2).name     = 'camk';
experiments(2).folder   = 'camk';
experiments(2).exclude  = {'ds220510a'}; % Removed at end of original script

% 3. PV
experiments(3).name     = 'pv';
experiments(3).folder   = 'PV';
experiments(3).exclude  = {'fr211206a', 'ds220728a'};

% --- Analysis Constants ---
SPECIAL_SUBJ_ID = 'ds220516a';        % Subject requiring custom channel map
CH_IDX_DEFAULT  = 6;                  % Standard start column for MUA
CH_IDX_SPECIAL  = [1, 4, 8, 10, 12, 14] + 5; % Custom columns for special subject

%% ========================================================================
%  2. MAIN PROCESSING LOOP
%  ========================================================================

for g = 1:length(experiments)
    exp_grp = experiments(g);
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Processing Group: %s\n', upper(exp_grp.name));
    fprintf('------------------------------------------------------------\n');
    
    % --- Define Paths ---
    dir_mua  = fullfile(PATH_MUA_ROOT, exp_grp.folder);
    dir_feat = fullfile(PATH_FEAT_ROOT, exp_grp.folder);
    
    % --- Get File Lists ---
    files_mua = dir(fullfile(dir_mua, '*.csv'));
    files_feat = dir(fullfile(dir_feat, '*.csv'));
    
    % --- Apply Exclusions ---
    if ~isempty(exp_grp.exclude)
        % Create a regex pattern for exclusion
        pattern = strjoin(exp_grp.exclude, '|');
        
        mask_mua = ~contains({files_mua.name}, exp_grp.exclude);
        files_mua = files_mua(mask_mua);
        
        mask_feat = ~contains({files_feat.name}, exp_grp.exclude);
        files_feat = files_feat(mask_feat);
        
        fprintf('Excluded subjects matching: %s\n', pattern);
    end
    
    % Sanity check: Ensure matched file counts
    if length(files_mua) ~= length(files_feat)
        warning('Mismatch in MUA and Feature file counts for %s', exp_grp.name);
    end
    
    % --- Load and Process Each Subject ---
    data_cells = cell(1, length(files_mua));
    sf_table_cells = cell(1, length(files_feat));
    
    for i = 1:length(files_mua)
        % 1. Load Tables
        tbl_mua  = readtable(fullfile(files_mua(i).folder, files_mua(i).name));
        tbl_feat = readtable(fullfile(files_feat(i).folder, files_feat(i).name));
        
        subj_id = unique(tbl_mua.id);
        if iscell(subj_id), subj_id = subj_id{1}; end
        
        % 2. Initialize Subject Data Table (Metadata)
        tbl_proc = tbl_mua(:, 1:5);
        
        % 3. Calculate Mean Firing Rate (FR)
        % Apply special channel selection if needed
        if strcmp(subj_id, SPECIAL_SUBJ_ID)
            cols_fr = CH_IDX_SPECIAL;
        else
            cols_fr = CH_IDX_DEFAULT:width(tbl_mua);
        end
        tbl_proc.FR_mean = mean(table2array(tbl_mua(:, cols_fr)), 2);
        
        % 4. Normalization (Z-score against Baseline)
        is_baseline = strcmp(tbl_feat.Condition, 'baseline');
        
        % A. Normalize Firing Rate
        base_mean_fr = mean(tbl_proc.FR_mean(is_baseline));
        base_std_fr  = std(tbl_proc.FR_mean(is_baseline));
        tbl_proc.FR_norm = (tbl_proc.FR_mean - base_mean_fr) / base_std_fr;
        
        % B. Normalize Spectral Slope Features (Cols 1:5)
        % Note: Slope is inverted per original logic
        tbl_feat.Slope = -tbl_feat.Slope; 
        
        feat_cols = 1:5;
        base_mean_feat = mean(tbl_feat{is_baseline, feat_cols}, 1);
        base_std_feat  = std(tbl_feat{is_baseline, feat_cols}, 0, 1);
        feats_norm = (tbl_feat{:, feat_cols} - base_mean_feat) ./ base_std_feat;
        
        % C. Normalize Power Spectrum Features (Cols 11:21)
        % Note: Slope fit is inverted per original logic
        tbl_feat.slope_fit = -tbl_feat.slope_fit;
        
        pow_cols = 11:21;
        base_mean_pow = mean(tbl_feat{is_baseline, pow_cols}, 1);
        base_std_pow  = std(tbl_feat{is_baseline, pow_cols}, 0, 1);
        pow_norm = (tbl_feat{:, pow_cols} - base_mean_pow) ./ base_std_pow;
        
        % 5. Consolidate Data
        % Append normalized arrays and the last column (original metadata/label)
        tbl_proc = [tbl_proc, array2table(feats_norm, 'VariableNames', tbl_feat.Properties.VariableNames(feat_cols)), ...
                              array2table(pow_norm, 'VariableNames', tbl_feat.Properties.VariableNames(pow_cols)), ...
                              tbl_feat(:, end)];
                          
        data_cells{i} = tbl_proc;
        sf_table_cells{i} = tbl_feat; % Keep raw spectral features for PSD extraction
    end
    
    % --- Aggregate Group Data ---
    DATA_all = vertcat(data_cells{:});
    DATA_sf_all = vertcat(sf_table_cells{:});
    
    % --- Extract Frequency and PSD Arrays ---
    % Columns 22 to end-1 contain [Freqs, PSDs] concatenated
    raw_psd_cols = 22:(width(DATA_sf_all)-1);
    n_pairs = length(raw_psd_cols) / 2;
    
    % Extract using efficient cell processing
    raw_data_mat = table2array(DATA_sf_all(:, raw_psd_cols));
    
    DATA_all.Frequencies = mat2cell(raw_data_mat(:, 1:n_pairs), ones(height(DATA_all),1));
    DATA_all.PSDs        = mat2cell(raw_data_mat(:, n_pairs+1:end), ones(height(DATA_all),1));
    
    % --- Save Results ---
    save_filename = sprintf('data_%s.mat', exp_grp.name);
    save_path = fullfile(PATH_OUT, save_filename);
    
    save(save_path, 'DATA_all');
    fprintf('Saved: %s\n', save_filename);
end

fprintf('\nData integration complete.\n');
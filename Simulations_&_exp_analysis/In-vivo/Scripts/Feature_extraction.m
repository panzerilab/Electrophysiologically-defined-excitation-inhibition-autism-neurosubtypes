% =========================================================================
% Title: Spectral Feature Extraction from chemogenetics manipulation
% experiments
%
% Author: Gabriele Mancini
% Date: 2025-12-10
%
% Description:
% This script processes experimental LFP and MUA neural data to extract periodic and
% aperiodic spectral features from LFP or EEG proxy signals. LFP/EEG signals
% are analyzed to quantify spectral properties (PSD slope, gamma power, Hurst exponent) alongside
% firing-rate statistics form MUA.
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

% Base Paths
ROOT_PATH = '../../In-silico';
DATA_ROOT = '../DATA/LFP';
SAVE_ROOT = '../results/features';

% Analysis Parameters
F_RANGE       = [30, 50];
F_RANGE_FOOOF = [20, 90];

% Add functions to path
addpath('../../functions');

%% ========================================================================
%  2. EXPERIMENT DEFINITIONS
%  ========================================================================
% Define structure array for experiments to avoid code duplication
% format: struct('name', 'folder_name_in_data', 'folder_name_in_results')

experiments(1).id = 'camk';
experiments(1).in_dir = fullfile(DATA_ROOT, 'camk');
experiments(1).out_dir = fullfile(SAVE_ROOT, 'camk');

experiments(2).id = 'hm4d';
experiments(2).in_dir = fullfile(DATA_ROOT, 'hm4d');
experiments(2).out_dir = fullfile(SAVE_ROOT, 'hm4d');

experiments(3).id = 'pv';
experiments(3).in_dir = fullfile(DATA_ROOT, 'PV');
experiments(3).out_dir = fullfile(SAVE_ROOT, 'PV');

%% ========================================================================
%  3. MAIN PROCESSING LOOP
%  ========================================================================

for g = 1:length(experiments)
    curr_exp = experiments(g);
    
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Processing Group: %s\n', upper(curr_exp.id));
    fprintf('------------------------------------------------------------\n');
    
    % 1. Create Output Directory
    if ~exist(curr_exp.out_dir, 'dir')
        mkdir(curr_exp.out_dir);
        fprintf('Created directory: %s\n', curr_exp.out_dir);
    end
    
    % 2. Get File List (Filtering hidden files and system dots)
    files = dir(curr_exp.in_dir);
    % Remove '.', '..', and files starting with '.' (like .DS_Store)
    valid_mask = ~ismember({files.name}, {'.', '..'}) & ~startsWith({files.name}, '.');
    file_list = files(valid_mask);
    
    % 3. Iterate through files
    for f = 1:length(file_list)
        fname = file_list(f).name;
        full_path = fullfile(curr_exp.in_dir, fname);
        
        fprintf('  [%d/%d] Loading: %s... ', f, length(file_list), fname);
        
        try
            % --- SPECIAL CASE LOGIC ---
            % Specific dataset requires channel removal
            if contains(fname, 'ds220516a')
                fprintf('(Special processing applied)\n');
                compute_features_remove_channels(full_path, curr_exp.out_dir, F_RANGE, F_RANGE_FOOOF);
            else
                compute_features(full_path, curr_exp.out_dir, F_RANGE, F_RANGE_FOOOF);
                fprintf('Done.\n');
            end
            
        catch ME
            fprintf('\n  Error processing %s: %s\n', fname, ME.message);
        end
    end
end

fprintf('\nAll processing complete.\n');
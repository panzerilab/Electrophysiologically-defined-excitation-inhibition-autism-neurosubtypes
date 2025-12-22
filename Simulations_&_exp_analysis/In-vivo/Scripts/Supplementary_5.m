% =========================================================================
% Title: Electrophysiology Regression & Longitudinal Modeling
%
% Author: Gabriele Mancini
% Date: 2025-12-10
%
% Description:
%   1. Loads LFP data (frequency corrected PSDs).
%   2. Performs subject-specific Z-score normalization against baseline 
%      (Log10 transformed).
%   3. Generates comparative PSD plots (Sham vs Treatment) for all groups.
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
%  1. CONFIGURATION & SETUP
%  ========================================================================

% --- Paths ---
ROOT_PATH = '../../In-silico';
SAVE_DIR  = '../Figures/Supplementary_5';
DATA_DIR  = '../results/tables'; % Assuming input files are here

% Create output directory
if ~exist(SAVE_DIR, 'dir')
    mkdir(SAVE_DIR);
    fprintf('Created directory: %s\n', SAVE_DIR);
end

addpath('../../functions');

% --- Constants ---
SMOOTH_WIN = 3;        % Moving average window for plotting
F_RANGE    = [0.1, 100]; % Hz

% --- Colors ---
C_SHAM   = [0.1059, 0.6196, 0.4667]; 
C_DREADD = [0.8510, 0.3725, 0.0078];

% --- Plotting Defaults ---
set(groot, ...
    'DefaultAxesFontSize', 10, ...
    'DefaultTextFontSize', 10, ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultTextFontName', 'Arial', ...
    'DefaultAxesLineWidth', 1.5, ...
    'DefaultAxesTickLabelInterpreter', 'none');

%% ========================================================================
%  2. DATA LOADING & PREPARATION
%  ========================================================================

% Define experiments struct to handle loading and labeling cleanly
exps(1).id = 'hm4di'; exps(2).id = 'camk'; exps(3).id = 'pv';
exps(1).file = 'data_hm4di.mat';
exps(2).file = 'data_camk.mat';
exps(3).file = 'data_pv.mat';

% Load and concatenate
data_cells = {};
for i = 1:length(exps)
    try
        % Load file (assuming variable structure matches previous scripts)
        tmp = load(fullfile(DATA_DIR, exps(i).file));
        
        % Check if variable is DATA_all or different
        if isfield(tmp, 'DATA_all')
            T = tmp.DATA_all;
        else
            % Fallback if variable name differs in _knee files
            fnames = fieldnames(tmp);
            T = tmp.(fnames{1});
        end
        
        % Add Labels
        T.exp_lab = repmat(string(exps(i).id), height(T), 1);
        
        data_cells{end+1} = T;
    catch ME
        warning('Could not load %s: %s', exps(i).file, ME.message);
    end
end

DATA = vertcat(data_cells{:});

% Apply Exclusion Filter (Specific PV Subject)
DATA = DATA(~(startsWith(DATA.id, 'fr') & DATA.exp_lab == "pv"), :);

%% ========================================================================
%  3. NORMALIZATION (Log10 -> Z-score vs Baseline)
%  ========================================================================
% 


fprintf('Computing baseline-normalized PSDs...\n');

% We process the whole table by iterating through unique subjects
DATA.PSD_norm = cell(height(DATA), 1);
subj_list = unique(DATA.id);

for s = 1:length(subj_list)
    curr_id = subj_list{s};
    
    % Get indices for this subject
    idx_subj = strcmp(DATA.id, curr_id);
    idx_base = idx_subj & strcmp(DATA.condition, 'baseline');
    
    % Compute Baseline Statistics
    if any(idx_base)
        % 1. Extract Baseline PSDs
        psd_base_raw = cell2mat(DATA.PSDs(idx_base));
        
        % 2. Log Transform
        psd_base_log = log10(psd_base_raw);
        
        % 3. Calculate Stats
        mu_base  = mean(psd_base_log, 1);
        sigma_base = std(psd_base_log, 0, 1);
        
        % 4. Apply Z-score to ALL trials for this subject
        row_indices = find(idx_subj);
        for r = row_indices'
            psd_curr = log10(DATA.PSDs{r});
            DATA.PSD_norm{r} = (psd_curr - mu_base) ./ sigma_base;
        end
    else
        warning('No baseline found for subject %s. Setting NaN.', curr_id);
        row_indices = find(idx_subj);
        for r = row_indices'
            DATA.PSD_norm{r} = NaN;
        end
    end
end

%% ========================================================================
%  4. PLOTTING
%  ========================================================================

fig = figure('Color', 'white', 'Units', 'centimeters', 'Position', [0, 0, 12, 16]);
t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Define specific plot order and titles
plot_order = {'hm4di', 'camk', 'pv'};
plot_titles = {'hSYN-hM4Di', 'CaMKII-hM3Dq', 'PV-hM4Di'};

for i = 1:3
    curr_exp = plot_order{i};
    nexttile; hold on;
    
    % --- Subset Data ---
    d_exp = DATA(DATA.exp_lab == curr_exp & ...
                 DATA.window_total <= 60 & ...
                 DATA.FR_norm <= 10 & ...
                 DATA.condition == "treatment", :);
             
    d_ctrl = d_exp(d_exp.grp == "ctrl", :);
    d_dread = d_exp(d_exp.grp == "exp", :);
    
    % --- Helper to Plot Group ---
    plot_group_psd(d_ctrl, C_SHAM,   SMOOTH_WIN, F_RANGE);
    plot_group_psd(d_dread, C_DREADD, SMOOTH_WIN, F_RANGE);
    
    % --- Formatting ---
    yline(0, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '--');
    set(gca, 'XScale', 'log', 'YScale', 'linear');
    xlim([0.2, 105]);
    
    title(plot_titles{i}, 'Interpreter', 'none', 'FontWeight', 'bold');
    ylabel('Norm. Power');
    grid off;
    
    if i == 1
        legend({'SHAM', 'DREADD'}, 'Location', 'best', 'Box', 'off');
    end
    if i == 3
        xlabel('Frequency (Hz)');
    else
        xticklabels([]);
    end
end

% Save
saveas(fig, fullfile(SAVE_DIR, 'PSDs_experiments.svg'));
fprintf('Figure saved to %s\n', fullfile(SAVE_DIR, 'PSDs_experiments.svg'));

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function plot_group_psd(tbl, color, win_smooth, f_lims)
    if isempty(tbl), return; end
    
    % Extract Freqs and PSDs
    freqs = double(tbl.Frequencies{1,:});
    psds  = cell2mat(tbl.PSD_norm);
    
    % Mask Frequencies
    mask = freqs >= f_lims(1) & freqs <= f_lims(2);
    f_plot = freqs(mask);
    p_plot = psds(:, mask);
    
    % Stats
    mu  = mean(p_plot, 1);
    sem = std(p_plot, 0, 1) ./ sqrt(size(p_plot, 1));
    
    % Smoothing
    mu_s  = movmean(mu, win_smooth);
    sem_s = movmean(sem, win_smooth);
    
    % Plot
    shadedErrorBar(f_plot, mu_s, sem_s, ...
        'lineprops', {'color', color, 'LineWidth', 1.5});
end
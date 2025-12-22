% =========================================================================
% Title: Electrophysiology Regression & Longitudinal Modeling
%
% Author: Gabriele Mancini
% Date: 2025-12-10
%
% Description:
%  Supplementary analysis pipeline focusing on single-feature regression
%  and group-specific prediction effects. The script loads electrophysiological
%  datasets across multiple genetic and chemogenetic manipulations, performs
%  artifact rejection and filtering, and constructs unified analysis tables.
%
%  Robust regression is applied separately to individual spectral features
%  (1/f slope and total power) to predict normalized firing rates, generating
%  feature-specific prediction maps. These predictions are then
%  evaluated within each experimental group (SHAM vs DREADD) and genetic
%  background (hSyn-hM4Di, CaMKII-hM3Dq, PV-hM4Di), producing subgroup-wise
%  correlation analyses and visualizations.
%
%  Finally, longitudinal mixed-effects modeling (LME) combined with
%  time-resolved Gramm plots assesses condition-by-group interactions
%  across recording windows for both predicted and raw features, providing
%  complementary temporal validation of the regression results. 
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
SAVE_DIR  = '../Figures/Supplementary_4';
DATA_DIR  = '../results/tables';

% Create output directory
if ~exist(SAVE_DIR, 'dir')
    mkdir(SAVE_DIR);
    fprintf('Created directory: %s\n', SAVE_DIR);
end

addpath('../../functions');

% --- Plotting Constants ---
FIG_WIDTH_SCATTER  = 10;
FIG_HEIGHT_SCATTER = 8;
FIG_WIDTH_PANEL    = 5;
FIG_HEIGHT_PANEL   = 17;

% Colors
C_HM4DI    = [117, 112, 179] / 255;
C_CAMKII   = [244, 166, 184] / 255;
C_PV       = [230, 171, 2] / 255;
C_SHAM     = [0.1059, 0.6196, 0.4667];
C_DREADD   = [0.8510, 0.3725, 0.0078];

% Set Default Figure Properties
set(groot, 'DefaultAxesFontSize', 10, 'DefaultTextFontSize', 10, ...
    'DefaultAxesFontName', 'Arial', 'DefaultAxesLineWidth', 1.5);

%% ========================================================================
%  2. DATA LOADING & PRE-PROCESSING
%  ========================================================================

try
    d_hm4di = load(fullfile(DATA_DIR, 'data_hm4di.mat'));
    d_camk  = load(fullfile(DATA_DIR, 'data_camK.mat'));
    d_pv    = load(fullfile(DATA_DIR, 'data_pv.mat'));
catch ME
    error('Failed to load data. Error: %s', ME.message);
end

% Assign Labels and Colors
process_data = @(S, lbl, col) addvars(S.DATA_all, ...
    repmat(string(lbl), height(S.DATA_all), 1), ...
    repmat(col, height(S.DATA_all), 1), ...
    'NewVariableNames', {'exp_lab', 'color'});

T_hm4di = process_data(d_hm4di, "hm4di", [0.4588, 0.4392, 0.7020]);
T_camk  = process_data(d_camk,  "camk",  [244, 166, 184] / 255);
T_pv    = process_data(d_pv,    "pv",    [230, 171, 2] / 255);

% Merge
DATA_all_all = vertcat(T_hm4di, T_camk, T_pv);

% Filters
% 1. Remove specific PV artifact
DATA_all_all = DATA_all_all(~(startsWith(DATA_all_all.id, 'fr') & DATA_all_all.exp_lab == "pv"), :);

% 2. Define Analysis Subset
filterData = @(data) data( ...
    data.window_total <= 60 & ...
    data.FR_norm <= 10 & ...
    data.condition == "treatment", :);

data_analysis = filterData(DATA_all_all);
y_target      = data_analysis.FR_norm;
group_colors  = data_analysis.color;

%% ========================================================================
%  3. PANELS E & F: SINGLE FEATURE REGRESSION
%  ========================================================================
% We loop through the two main single features: Slope (Panel E) and Tot_Pow (Panel F)

reg_tasks = {
    'Slope',   'panel_E.svg', 'pred_slope';
    'tot_pow', 'panel_F.svg', 'pred_totp'
};

fprintf('Generating Panels E and F...\n');

for k = 1:size(reg_tasks, 1)
    feat_name = reg_tasks{k, 1};
    save_name = reg_tasks{k, 2};
    pred_col  = reg_tasks{k, 3};
    
    % Prepare Figure
    fig = figure('Units', 'centimeters', 'Position', [0, 0, FIG_WIDTH_SCATTER, FIG_HEIGHT_SCATTER], 'Color', 'w');
    hold on;
    
    X = data_analysis{:, feat_name};
    
    % Robust Regression
    [B_c, ~, y_pred, y_fit, ~, ~, r, ~] = performRegression_rob(X, y_target, 100);
    
    % Store prediction in main table for later panels
    % Note: We calculate prediction for the WHOLE dataset, not just the filtered subset
    DATA_all_all.(pred_col) = B_c(1) + B_c(2) .* DATA_all_all.(feat_name);
    
    % Plot
    scatter(y_target, y_pred, 10, group_colors, 'filled');
    plot(y_target, y_fit, 'k', 'LineWidth', 2);
    
    xlabel('Actual FR (z-scored)');
    ylabel('Predicted FR (z-scored)');
    
    % Annotations (Dynamic positioning based on data range could be added here)
    if strcmp(feat_name, 'Slope')
        text(0, -1, sprintf('Slope: r = %.2f ***', r), 'FontWeight', 'bold');
        y_txt = -2.4; 
    else
        text(1, -2, sprintf('Total Power: r = %.2f ***', r), 'FontWeight', 'bold');
        y_txt = -4;
    end
    
    text(1, y_txt,       'hSYN-hM4Di',   'Color', C_HM4DI,  'FontWeight', 'bold');
    text(1, y_txt - 0.6, 'CaMKII-hM3Dq', 'Color', C_CAMKII, 'FontWeight', 'bold');
    text(1, y_txt - 1.2, 'PV-hM4Di',     'Color', C_PV,     'FontWeight', 'bold');
    
    saveas(fig, fullfile(SAVE_DIR, save_name));
    close(fig);

end

%% ========================================================================
%  4. PANELS C & D: PREDICTION BY GROUP (SUBPLOTS)
%  ========================================================================
% Loop through Slope Predictions (Panel C) and Tot_Pow Predictions (Panel D)

panel_tasks = {
    'pred_slope', 'panel_C.svg';
    'pred_totp',  'panel_D.svg'
};

exp_labs_list = {"hm4di", "camk", "pv"};

fprintf('Generating Panels C and D...\n');

for k = 1:size(panel_tasks, 1)
    pred_var = panel_tasks{k, 1};
    save_name = panel_tasks{k, 2};
    
    fig = figure('Units', 'centimeters', 'Position', [0, 0, FIG_WIDTH_PANEL, FIG_HEIGHT_PANEL], 'Color', 'w');
    
    for i = 1:3
        exp_name = exp_labs_list{i};
        
        % Subplot: Top=HM4Di, Mid=CaMK, Bot=PV
        subplot(3, 1, i); hold on;
        
        % --- Subset Data ---
        mask_base = DATA_all_all.window_total <= 60 & ...
                    DATA_all_all.FR_norm <= 10 & ...
                    DATA_all_all.condition == "treatment" & ...
                    DATA_all_all.exp_lab == exp_name;
        
        d_ctrl = DATA_all_all(mask_base & DATA_all_all.grp == "ctrl", :);
        d_exp  = DATA_all_all(mask_base & DATA_all_all.grp == "exp", :);
        
        % --- Prepare Vectors ---
        FR_ctrl = d_ctrl.FR_norm;   P_ctrl = d_ctrl.(pred_var);
        FR_exp  = d_exp.FR_norm;    P_exp  = d_exp.(pred_var);
        
        % --- Correlations ---
        [R_c, P_val_c] = corr(FR_ctrl, P_ctrl, 'Type', 'Spearman');
        [R_e, P_val_e] = corr(FR_exp, P_exp, 'Type', 'Spearman');
        
        % --- Plotting ---
        % Sham
        scatter(FR_ctrl, P_ctrl, 5, 'o', 'MarkerEdgeColor', C_SHAM, 'MarkerFaceColor', C_SHAM);
        p1 = polyfit(FR_ctrl, P_ctrl, 1);
        [x1_sort, idx1] = sort(FR_ctrl); % Sort for clean line
        plot(x1_sort, polyval(p1, x1_sort), 'Color', C_SHAM, 'LineWidth', 1.5);
        
        % DREADD
        scatter(FR_exp, P_exp, 5, 'o', 'MarkerEdgeColor', C_DREADD, 'MarkerFaceColor', C_DREADD);
        p2 = polyfit(FR_exp, P_exp, 1);
        [x2_sort, idx2] = sort(FR_exp);
        plot(x2_sort, polyval(p2, x2_sort), 'Color', C_DREADD, 'LineWidth', 1.5);
        
        % Formatting
        if i == 3, xlabel('Actual FR (z-scored)'); else, xticklabels([]); end
        ylabel('Predicted FR (z-scored)');
        
        % Stats Text
        y_min = min(P_exp);
        text(-4, 0.95*y_min, sprintf('SHAM: r=%.2f p=%.3f', R_c, P_val_c), 'Color', C_SHAM, 'FontSize', 8, 'FontWeight', 'bold');
        text(-4, 0.55*y_min, sprintf('DREADD: r=%.2f p=%.3f', R_e, P_val_e), 'Color', C_DREADD, 'FontSize', 8, 'FontWeight', 'bold');
    end
    
    saveas(fig, fullfile(SAVE_DIR, save_name));
    close(fig);
end

%% ========================================================================
%  5. LONGITUDINAL ANALYSIS (LME + GRAMM)
%  ========================================================================
% Performs LME and plots time-series for: 
% 1. Predicted Slope, 2. Predicted Tot_Pow, 3. Raw Slope

long_tasks = {
    'pred_slope', 'Predicted FR', 'Predictions_SLOPE';
    'pred_totp',  'Predicted FR', 'Predictions_TOT';
    'Slope',      'Slope',        'Slope'
};

fprintf('Generating Longitudinal plots...\n');

for t = 1:size(long_tasks, 1)
    
    curr_var   = long_tasks{t, 1};
    curr_label = long_tasks{t, 2};
    curr_tag   = long_tasks{t, 3};
    
    for i = 1:3
        exp_name = exp_labs_list{i};
        
        % Filter Data
        d_sub = DATA_all_all(DATA_all_all.window_total <= 60 & ...
                             DATA_all_all.FR_norm <= 10 & ...
                             DATA_all_all.exp_lab == exp_name, :);
                         
        % Standardize
        d_sub.time = d_sub.window_total;
        d_sub.condition(strcmp(d_sub.condition, 'cno')) = {'treatment'};
        
        % --- LME Model ---
        formula = sprintf('%s ~ grp*condition + (1|id)', curr_var);
        try
            mod = fitlme(d_sub, formula);
            anova_res = anova(mod);
            p_val = anova_res.pValue(4);
            fprintf('[%s] %s Interaction P-value: %.4f\n', exp_name, curr_var, p_val);
        catch
            fprintf('LME failed for %s in %s\n', curr_var, exp_name);
        end
        
        % --- Gramm Plot ---
        clear g;
        g = gramm('x', d_sub.window_total, 'y', d_sub.(curr_var), 'color', d_sub.grp);
        g.stat_smooth();
        g.geom_vline('xintercept', 15);
        g.geom_vline('xintercept', 35);
        g.set_names('x', 'Time Window', 'y', curr_label);
        
        % Create Figure
        fig = figure('Units', 'centimeters', 'Position', [0, 0, 6, 5], 'Color', 'w');
        g.draw();
        
        % Save
        fname = sprintf('%s_vs_time_%s.svg', curr_tag, exp_name);
        saveas(fig, fullfile(SAVE_DIR, fname));
        close(fig);
    end
end

fprintf('Supplementary Analysis Complete.\n');
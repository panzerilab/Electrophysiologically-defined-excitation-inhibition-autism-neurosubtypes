% =========================================================================
% Title: Electrophysiology Regression & Longitudinal Modeling
%
% Author: Gabriele Mancini
% Date: 2025-12-10
%
% Description:
%  End-to-end analysis pipeline for electrophysiology-based regression and
%  longitudinal modeling. The script loads multi-unit activity (MUA) and
%  LFP-derived spectral features across experimental conditions and genetic
%  manipulations, performs artifact rejection and data filtering, applies
%  normalization relative to baseline periods, and aggregates data into
%  subject- and experiment-specific tables.
%
%  The pipeline implements k-fold cross-validated linear regression and
%  robust regression to compare multiple feature sets (e.g., 1/f slope,
%  total power, H, gamma, and H + gamma) in predicting normalized firing
%  rates. Model performance is quantified via correlation metrics, paired
%  statistical tests, and extensive bootstrap procedures.
%
%  Additional analyses include simulated-model validation, group-wise
%  predicted-versus-actual comparisons, and longitudinal mixed-effects
%  modeling (LME) to assess condition-by-group interactions over time,
%  with effect sizes reported. All results are summarized in publication-
%  ready figures and exported tables to support reproducibility.
%
% Associated manuscript:
% Bertelsen N. et al.,
% "Electrophysiologically-defined excitation–inhibition autism neurosubtypes",
% Nature Neuroscience (under review).
%
% License: MIT License
% Contact: gabriele.mancini@oiit.it
% =========================================================================


close all; clc;
% ---------------------------
% CONFIG
% ---------------------------
ROOT_PATH = '../../In-silico';
SAVE_DIR  = '../Figures/Figure 2';
DATA_DIR  = '../results/tables';

% Add custom functions path
p = genpath('../../functions');
p = strrep(p, [pathsep '../../functions/.git'], '');
addpath(p);

RANDOM_SEED = 234;
K_FOLDS = 10;
N_BOOT_BOOTSTRAP = 1000; % for bootstrap z-tests
N_BOOT_ROBUST = 100;     % used by robust regression where applicable

% plotting constants (centimeters)
FIG_WIDTH_STD  = 15;
FIG_HEIGHT_STD = 7.5;

% color palette
C_PASTEL_RED   = [0.4588, 0.4392, 0.7020];
C_PASTEL_GREEN = [244, 166, 184] / 255;
C_PASTEL_BLUE  = [230, 171, 2] / 255;
C_SHAM         = [0.1059, 0.6196, 0.4667];
C_DREADD       = [0.8510, 0.3725, 0.0078];
C_HM4DI        = [117, 112, 179] / 255;
C_CAMKII       = [244, 166, 184] / 255;
C_PV           = [230, 171, 2] / 255;

% Ensure SAVE_DIR exists
if ~exist(SAVE_DIR, 'dir')
    mkdir(SAVE_DIR);
    fprintf('Created directory: %s\n', SAVE_DIR);
end

% set rng
rng(RANDOM_SEED);

% set default figure properties
set(groot, ...
    'DefaultAxesFontSize', 10, ...
    'DefaultTextFontSize', 10, ...
    'DefaultLegendFontSize', 10, ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultTextFontName', 'Arial', ...
    'DefaultLegendFontName', 'Arial', ...
    'DefaultAxesLineWidth', 1.5, ...
    'DefaultAxesTickLabelInterpreter', 'tex');

% ---------------------------
% LOAD & PROCESS DATA
% ---------------------------
try
    d_hdm4di = load(fullfile(DATA_DIR, 'data_hm4di.mat'));
    d_camk   = load(fullfile(DATA_DIR, 'data_camK.mat'));
    d_pv     = load(fullfile(DATA_DIR, 'data_pv.mat'));
catch ME
    error('Failed to load datasets from %s. Error: %s', DATA_DIR, ME.message);
end

% helper: add exp label and color variable
process_table = @(T, lab, col) addvars(T.DATA_all, ...
    repmat(string(lab), height(T.DATA_all), 1), ...
    repmat(col, height(T.DATA_all), 1), ...
    'NewVariableNames', {'exp_lab', 'color'});

T_hdm4di = process_table(d_hdm4di, "hm4di", C_PASTEL_RED);
T_camk   = process_table(d_camk,   "camk",  C_PASTEL_GREEN);
T_pv     = process_table(d_pv,     "pv",    C_PASTEL_BLUE);

DATA_all_all = vertcat(T_hdm4di, T_camk, T_pv);

% Remove artifact: specific PV/FR id pattern
if any(contains(DATA_all_all.Properties.VariableNames, 'id'))
    DATA_all_all = DATA_all_all(~(startsWith(DATA_all_all.id, 'fr') & DATA_all_all.exp_lab == "pv"), :);
end

% filter function for analysis
get_valid_data = @(df) df( ...
    df.window_total <= 60 & ...
    df.FR_norm <= 10 & ...
    df.condition == "treatment", :);

data_analysis  = get_valid_data(DATA_all_all);

% ---------------------------
% DEFINE FEATURES & TARGET
% ---------------------------
features = { "Slope", "tot_pow", ["H", "gamma"], "H", "gamma" }; % cell array; 3rd is combined
n_features = numel(features);
y_target = data_analysis.FR_norm;

% convert features into X matrices on demand in modular parts
% ---------------------------
% PANEL C: regression analysis & model comparison
% ---------------------------
fprintf('Running regression analyses (Panel C)...\n');

% initialize containers
real_corr = zeros(1, n_features);
real_corr_sem = zeros(1, n_features);
all_r_boots = nan(K_FOLDS, n_features + 1); % +1 for sim model appended later
best_model_weights = []; % will be set for the combined H+gamma model

for i = 1:n_features
    feat_set = features{i};
    X_mat = data_analysis{:, feat_set}; % works for single or multi-column sets
    try
        [B_weights, r_boot] = performRegression_kfold(X_mat, y_target, K_FOLDS, RANDOM_SEED);
    catch ME
        warning('performRegression_kfold failed for feature set %s: %s', strjoin(string(feat_set), '+'), ME.message);
        B_weights = zeros(size(X_mat,2)+1,1); % fallback: zeros
        r_boot = zeros(K_FOLDS,1);
    end

    % mean & sem across folds
    real_corr(i) = mean(r_boot);
    real_corr_sem(i) = std(r_boot, 0, 1) / sqrt(K_FOLDS);
    all_r_boots(1:K_FOLDS, i) = r_boot(:);

    % Save best model weights for combined H+gamma (3rd entry originally)
    if isequal(string(feat_set), string(features{3}))
        B_c_mean = mean(B_weights, 1);    % if B_weights is (folds x params)
        if isrow(B_c_mean), B_c_mean = B_c_mean(:)'; end
        best_model_weights = B_c_mean;    % store for Panel D prediction
        % Persist a lightweight MAT file for reproducibility
        try
            save(fullfile(SAVE_DIR, 'Weights_all_rob.mat'), 'B_c_mean');
        catch
            warning('Could not save Weights_all_rob.mat');
        end
    end
end

% B. K-Fold Validation Simulation (Simulated Model)
feat_sim = features{3}; 
X_sim = data_analysis{:, feat_sim};
y_sim = y_target(:);

cv = cvpartition(length(y_sim), 'KFold', K_FOLDS);
sim_corrs = zeros(1, K_FOLDS);

for k = 1:K_FOLDS
    tr_idx = training(cv, k);
    te_idx = test(cv, k);
    
    B_fold = (X_sim(tr_idx,:)' * X_sim(tr_idx,:)) \ (X_sim(tr_idx,:)' * y_sim(tr_idx));
    y_p_fold = X_sim(te_idx,:) * B_fold;
    sim_corrs(k) = corr(y_sim(te_idx), y_p_fold, 'Type', 'Pearson');
end

mean_R_sim = mean(sim_corrs);
sem_R_sim = std(sim_corrs) / sqrt(K_FOLDS);
real_corr = [real_corr, mean_R_sim];
real_corr_sem = [real_corr_sem, sem_R_sim];
all_r_boots(:, n_features+1) = sim_corrs(:);

% Create labels & variable names for tables & plotting
feature_labels = [features, {"sim"}];
lbls_display = cellfun(@(x) strjoin(string(x), ' + '), feature_labels, 'UniformOutput', false);
lbls_display = string(lbls_display);
lbls_vars = arrayfun(@(s) matlab.lang.makeValidName(replace(s, ' + ', '_plus_')), lbls_display);

% Plot comparison bar (Panel C simplified)
fig_c = figure('Units','centimeters','Position',[0,0,FIG_WIDTH_STD,FIG_HEIGHT_STD],'Color','w');
subplot(1,3,1);
% If best_model_weights exists, plot it
if ~isempty(best_model_weights)
    B_c_mean = best_model_weights;
    % If B_c_mean includes intercept as first element, plot only weights 2:end
    wvals = B_c_mean(2:end);
    wstd = zeros(size(wvals));
    if size(B_c_mean,2) > 1
        % Attempt to compute std if data available (if B_weights was folds x params)
        % Here we don't have B_weights; keep wstd zeros to avoid crash
        wstd = zeros(size(wvals));
    end
    bar(wvals, 'FaceColor', [0.7 0.7 0.7]);
    hold on;
    errorbar(1:length(wvals), wvals, wstd, 'k.', 'LineWidth', 0.8, 'CapSize', 10);
    xticks(1:length(wvals));
    xticklabels({'H', '\gamma'});
    ylabel('Regression Weight');
    title('Model Weights (H + gamma)');
    box off;
else
    text(0.5,0.5,'No weights available','HorizontalAlignment','center');
end

% Panel: model comparison bar
subplot(1,3,[2 3]);
hold on;
bar(real_corr, 'FaceColor', 'flat', 'EdgeColor', 'none');
errorbar(1:numel(real_corr), real_corr, real_corr_sem, 'k.', 'LineWidth', 1);
xticks(1:numel(real_corr));
xticklabels(lbls_display);
xtickangle(45);
ylabel('Regression Performance r');
title('Model Comparison');
ylim([-0.1, 0.9]);
box off;
saveas(fig_c, fullfile(SAVE_DIR, 'panel_JK.svg'));

% ---------------------------
% Pairwise t-tests matrix across models (using per-fold correlations)
% ---------------------------
n_models = size(all_r_boots, 2);
p_matrix = nan(n_models);
for a = 1:n_models
    for b = a+1:n_models
        try
            [~, p_val] = ttest(all_r_boots(:,a), all_r_boots(:,b));
        catch
            p_val = NaN;
        end
        p_matrix(a,b) = p_val;
        p_matrix(b,a) = p_val;
    end
end
p_table = array2table(round(p_matrix,3), 'VariableNames', lbls_vars, 'RowNames', lbls_display);

% Save pairwise pvalues CSV
writetable(p_table, fullfile(SAVE_DIR, 'model_pvalues.csv'), 'WriteRowNames', true);

% ---------------------------
% PANEL I: Best model scatter & robust regression plot
% ---------------------------
fprintf('Panel I (scatter) and robust regression...\n');
fig_i = figure('Units','centimeters','Position',[0,0,10,8],'Color','w');
X_best = data_analysis{:, features{3}};
try
    [~, ~, y_pred_rob, y_fit_rob, ~, ~, r_rob, ~] = performRegression_rob(X_best, y_target, N_BOOT_ROBUST);
catch ME
    warning('performRegression_rob failed: %s', ME.message);
    y_pred_rob = zeros(size(y_target));
    y_fit_rob = zeros(size(y_target));
    r_rob = NaN;
end
scatter(y_target, y_pred_rob, 10, data_analysis.color, 'filled');
hold on;
plot(y_target, y_fit_rob, 'k', 'LineWidth', 2);
xlabel('Actual FR (z-scored)');
ylabel('Predicted FR (z-scored)');
text(0.01, 0.98, sprintf('Best Model (H + \\gamma): r = %.2f', r_rob), 'Units','normalized', 'VerticalAlignment','top', 'FontWeight','bold');
saveas(fig_i, fullfile(SAVE_DIR, 'panel_I.svg'));

% ---------------------------
% PANEL D: Predicted vs Actual by group
% ---------------------------
fprintf('Panel D: predicted vs actual by group...\n');
if isempty(best_model_weights)
    warning('best_model_weights missing. Predictions will be zero.');
    % create zero weights with intercept + H + gamma placeholders
    best_model_weights = [0, 0, 0];
end

% Ensure required columns exist in DATA_all_all; create pred variable
if ~ismember('pred', DATA_all_all.Properties.VariableNames)
    % assume ordering: [intercept, H, gamma]
    DATA_all_all.pred = best_model_weights(1) + ...
        best_model_weights(2) .* DATA_all_all.H + ...
        best_model_weights(3) .* DATA_all_all.gamma;
end

fig_d = figure('Units','centimeters','Position',[0,0,10,17],'Color','w');
exp_labs_list = {"hm4di", "camk", "pv"};
panelD_stats = cell2table(cell(0,7), ...
    'VariableNames', {'Experiment','R_ctrl','P_ctrl','R_exp','P_exp','R_all_z','P_all_z'});


for i = 1:3
    exp_name = exp_labs_list{i};
    mask_base = DATA_all_all.window_total <= 60 & ...
                DATA_all_all.FR_norm <= 10 & ...
                DATA_all_all.condition == "treatment" & ...
                DATA_all_all.exp_lab == exp_name;

    d_ctrl = DATA_all_all(mask_base & DATA_all_all.grp == "ctrl", :);
    d_exp  = DATA_all_all(mask_base & DATA_all_all.grp == "exp", :);

    FR_ctrl = d_ctrl.FR_norm;   pred_ctrl = d_ctrl.pred;
    FR_exp  = d_exp.FR_norm;    pred_exp  = d_exp.pred;

    FR_all      = [FR_ctrl; FR_exp];
    pred_all    = [pred_ctrl; pred_exp];

    FR_all_cen   = FR_all - mean(FR_all);
    pred_all_cen = pred_all - mean(pred_all);

    % Correlations (Spearman)
    try
        [R_c, P_c] = corr(FR_ctrl, pred_ctrl, 'Type', 'Spearman');
    catch
        R_c = NaN; P_c = NaN;
    end
    try
        [R_e, P_e] = corr(FR_exp, pred_exp, 'Type', 'Spearman');
    catch
        R_e = NaN; P_e = NaN;
    end
    try
        [R_cen, P_cen] = corr(FR_all_cen, pred_all_cen, 'Type', 'Spearman');
    catch
        R_cen = NaN; P_cen = NaN;
    end

    panelD_stats(end+1,:) = {exp_name, R_c, P_c, R_e, P_e, R_cen, P_cen};

    % Plot left: raw correlations (Ctrl, Exp)
    subplot(3,2,2*(i-1)+1); hold on;
    scatter(FR_ctrl, pred_ctrl, 5, 'o', 'MarkerEdgeColor', C_SHAM, 'MarkerFaceColor', C_SHAM);
    if ~isempty(FR_ctrl)
        p1 = polyfit(FR_ctrl, pred_ctrl, 1);
        plot(FR_ctrl, polyval(p1, FR_ctrl), 'Color', C_SHAM, 'LineWidth', 1.5);
    end

    scatter(FR_exp, pred_exp, 5, 'o', 'MarkerEdgeColor', C_DREADD, 'MarkerFaceColor', C_DREADD);
    if ~isempty(FR_exp)
        p2 = polyfit(FR_exp, pred_exp, 1);
        plot(FR_exp, polyval(p2, FR_exp), 'Color', C_DREADD, 'LineWidth', 1.5);
    end

    p_all = polyfit(FR_all, pred_all, 1);
    plot(FR_all, polyval(p_all, FR_all), 'k', 'LineWidth', 1.5);

    ylabel('Predicted FR');
    if i == 3, xlabel('Actual FR'); end

    % Annotations
    y_max = max([FR_ctrl; FR_exp; 1]);
    text(-4, 0.95*y_max, sprintf('SHAM: r=%.2f p=%.3f', R_c, P_c), 'Color', C_SHAM, 'FontSize', 8);
    text(-4, 0.55*y_max, sprintf('DREADD: r=%.2f p=%.3f', R_e, P_e), 'Color', C_DREADD, 'FontSize', 8);

    % Right: centered correlations
    subplot(3,2,2*i); hold on;
    scatter(FR_all_cen, pred_all_cen, 10, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', C_SHAM);
    p_cen = polyfit(FR_all_cen, pred_all_cen, 1);
    plot(FR_all_cen, polyval(p_cen, FR_all_cen), 'k', 'LineWidth', 1.5);
    ylabel('Pred FR (z)');
    if i == 3, xlabel('Actual FR (z)'); end
    text(-4, mean(FR_all_cen)+3, sprintf('All: r=%.2f p=%.3f', R_cen, P_cen), 'Color', 'k', 'FontSize', 8);
end
saveas(fig_d, fullfile(SAVE_DIR, 'panel_D.svg'));

% Convert panelD_stats to table & save
writetable(panelD_stats, fullfile(SAVE_DIR,'panelD_stats.csv'));

%% ========================================================================
%  6. LONGITUDINAL ANALYSIS (LME & Gramm Plots) with Effect Size Table
% ========================================================================

fprintf('Generating Longitudinal Plots and LME statistics...\n');

metrics = {'pred', 'H', 'gamma_pow', 'FR_norm'};
metric_labels = {'Predicted FR', 'H', 'Gamma Pow.', 'FR norm'};

YLIMS = cell(4,3);
YLIMS{1,1} = [-5 5];   YLIMS{1,2} = [-2 8];   YLIMS{1,3} = [-4 10];
YLIMS{4,1} = [-5 5];   YLIMS{4,2} = [-2 8];   YLIMS{4,3} = [-4 10];

% Initialize results
LME_results = {}; % Experiment | Metric | F | DF1 | DF2 | p | eta_p2

for exp_idx = 1:3
    exp_name = exp_labs_list{exp_idx};
    
    % Filter data
    d_sub = DATA_all_all(DATA_all_all.window_total <= 60 & ...
                         DATA_all_all.FR_norm <= 10 & ...
                         DATA_all_all.exp_lab == exp_name, :);
    
    % Standardize condition
    if iscell(d_sub.condition)
        d_sub.condition = string(d_sub.condition);
    end
    d_sub.condition(strcmp(d_sub.condition, 'cno')) = "treatment";
    d_sub.time = d_sub.window_total;
    
    for m_idx = 1:length(metrics)
        met = metrics{m_idx};
        formula = sprintf('%s ~ grp*condition + (1|id)', met);
        
        try
            % Fit LME
            lme_mod = fitlme(d_sub, formula);
            
            % ANOVA table (Type III marginal tests)
            anova_res = anova(lme_mod, 'DFMethod', 'Residual');
            
            % Find interaction term
            interaction_row = find(strcmp(anova_res.Term, 'condition:grp'));
            if ~isempty(interaction_row)
                Fval = anova_res.FStat(interaction_row);
                DF1  = anova_res.DF1(interaction_row);
                DF2  = anova_res.DF2(interaction_row);
                p_val = anova_res.pValue(interaction_row);
                eta_p2 = (Fval * DF1) / (Fval * DF1 + DF2);
            else
                warning('Interaction term not found in ANOVA table.');
                Fval = NaN; DF1 = NaN; DF2 = NaN; p_val = NaN; eta_p2 = NaN;
            end
            
            % Store results
            LME_results(end+1,:) = {exp_name, met, Fval, DF1, DF2, p_val, eta_p2};
            
            
        catch ME
            warning('LME failed for %s - %s: %s', exp_name, met, ME.message);
            LME_results(end+1,:) = {exp_name, met, NaN, NaN, NaN, NaN, NaN};
        end
        
        % --- Gramm Plot ---
        try
            clear g;
            g = gramm('x', d_sub.window_total, 'y', d_sub.(met), 'color', d_sub.grp);
            g.stat_smooth();
            g.geom_vline('xintercept', 15);
            g.geom_vline('xintercept', 35);
            g.set_names('x', 'Time Window', 'y', metric_labels{m_idx});
            if ~isempty(YLIMS{m_idx, exp_idx})
                g.axe_property('YLim', YLIMS{m_idx, exp_idx});
            end
            fig_g = figure('Units','centimeters','Position',[0,0,6,5],'Color','w');
            g.draw();
            fname = sprintf('%s_%s.svg', met, exp_name);
            if strcmp(met,'gamma_pow')
                fname = sprintf('gamma_%s.svg', exp_name);
            end
            saveas(fig_g, fullfile(SAVE_DIR, fname));
            close(fig_g);
        catch ME
            warning('Gramm plotting failed for %s - %s: %s', exp_name, met, ME.message);
        end
    end
end

% Convert to table
LME_tbl = cell2table(LME_results, ...
    'VariableNames', {'Experiment','Metric','F','DF1','DF2','p_interaction','eta_p2'});

% Save table
writetable(LME_tbl, fullfile(SAVE_DIR, 'stats_LME_interaction_effectsizes.csv'));

% Display table in MATLAB terminal
disp(LME_tbl);


%% ---------------------------
% BOOTSTRAPPED COMPARISONS (robust r bootstrap & pairwise z-tests)
% ---------------------------
fprintf('Running bootstrap comparisons...\n');
n_boot = N_BOOT_BOOTSTRAP;
r_stats_mean = zeros(n_features, 1);
r_stats_boot = zeros(n_boot, n_features);

for i = 1:n_features
    feat_i = features{i};
    X_f = data_analysis{:, feat_i};
    y_f = data_analysis.FR_norm;
    r_vals = zeros(n_boot,1);

    for b = 1:n_boot
        idx = randi(size(X_f,1), size(X_f,1), 1);
        try
            [~, ~, ~, ~, ~, ~, r_curr, ~] = performRegression_rob(X_f(idx,:), y_f(idx), 1);
        catch
            r_curr = NaN;
        end
        r_vals(b) = r_curr;
    end
    r_stats_boot(:, i) = r_vals;
    r_stats_mean(i) = nanmean(r_vals);
end

% pairwise z-tests
comp_results = {};
for i = 1:n_features
    for j = i+1:n_features
        m1 = r_stats_mean(i);
        m2 = r_stats_mean(j);
        se1 = nanstd(r_stats_boot(:, i));
        se2 = nanstd(r_stats_boot(:, j));
        se_pooled = sqrt(se1^2 + se2^2);
        if se_pooled == 0 || isnan(se_pooled)
            z_stat = NaN;
            p_val_z = NaN;
        else
            z_stat = (m1 - m2) / se_pooled;
            p_val_z = 2 * (1 - normcdf(abs(z_stat)));
        end
        comp_results(end+1, :) = { strjoin(features{i}, '+'), strjoin(features{j}, '+'), m1, m2, z_stat, p_val_z };
    end
end

results_table = cell2table(comp_results, 'VariableNames', {'Feature1','Feature2','Mean_r1','Mean_r2','Z_stat','p_value'});
writetable(results_table, fullfile(SAVE_DIR, 'bootstrap_comparisons.csv'));

% ---------------------------
% COLLECT & EXPORT ALL STATS (master tables)
% ---------------------------
fprintf('Exporting summary tables...\n');
% 1) Regression performance (numeric & formatted)
r_mean = real_corr(:);
r_sem = real_corr_sem(:);
feat_names = cellstr(lbls_display(:));

T_numeric = table(r_mean, r_sem, 'VariableNames', {'r_mean','r_sem'}, 'RowNames', feat_names);
writetable(T_numeric, fullfile(SAVE_DIR, 'regression_performance_numeric.csv'), 'WriteRowNames', true);

fmt_strings = arrayfun(@(m,s) sprintf('%.3f +/- %.3f', m, s), r_mean, r_sem, 'UniformOutput', false);
T_formatted = table(fmt_strings, 'VariableNames', {'r_mean_sem'}, 'RowNames', feat_names);
writetable(T_formatted, fullfile(SAVE_DIR, 'regression_performance_formatted.csv'), 'WriteRowNames', true);

% 2) Pairwise t-tests (already have p_table)
writetable(p_table, fullfile(SAVE_DIR, 'stats_pairwise_ttests.csv'), 'WriteRowNames', true);

% 3) Bootstrapped z-tests (results_table)
writetable(results_table, fullfile(SAVE_DIR, 'stats_bootstrap_ztests.csv'));

% 4) Panel D correlations (already saved earlier as panelD_tbl)
% 5) LME interactions table (already saved)

% 6) Save MATLAB .mat snapshot with key variables
try
    save(fullfile(SAVE_DIR, 'analysis_snapshot.mat'), 'real_corr', 'real_corr_sem', 'all_r_boots', ...
        'lbls_display', 'lbls_vars', 'results_table', 'best_model_weights', 'features');
catch
    warning('Could not save analysis snapshot.');
end

fprintf('Analysis complete. All figures and tables are in: %s\n', SAVE_DIR);


%% ------------------------------------------------------------
% DISPLAY ALL TABLES IN MATLAB TERMINAL
% ------------------------------------------------------------

fprintf('\n===== PANEL B STATISTICS =====\n');
disp(LME_tbl);

fprintf('\n===== PANEL C STATISTICS =====\n');
disp(p_table);

fprintf('\n===== PANEL D STATISTICS =====\n');
disp(panelD_stats);

fprintf('\n===== PANEL E STATISTICS =====\n');
disp(results_table);

% Ensure column vectors
r_mean = real_corr(:);
r_sem  = real_corr_sem(:);
feat_names = cellstr(lbls_display(:)); % cell array for row names

% 1) Numeric table (easy to use for downstream stats)
T_numeric = table(r_mean, r_sem, ...
    'VariableNames', {'r_mean', 'r_sem'}, ...
    'RowNames', feat_names);

% Display numeric table
disp('Numeric regression performance (mean, sem):');
disp(T_numeric);

% Save numeric table (RowNames included)
writetable(T_numeric, fullfile(SAVE_DIR, 'regression_performance_numeric.csv'), ...
    'WriteRowNames', true);

% 3) (Optional) also save as MATLAB .mat for later use
save(fullfile(SAVE_DIR, 'regression_performance_tables.mat'), 'T_numeric');


fprintf('\nAll statistical tables have been printed to the MATLAB terminal.\n');

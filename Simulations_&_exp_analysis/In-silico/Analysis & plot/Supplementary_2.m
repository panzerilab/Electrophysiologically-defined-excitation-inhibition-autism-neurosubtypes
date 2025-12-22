% =========================================================================
% Title: Statistical Analysis and Plotting of Spectral Features
%
% Author: Gabriele Mancini
% Date: 2025-12-10
%
% Description:
% This script loads precomputed spectral and firing-rate features extracted
% from simulated LFP or EEG proxy signals. It performs statistical analyses
% (correlation, cross-validated regression) and generates the panels for
% Supplementary Figure 1, including predicted-versus-actual firing-rate
% plots and representative power spectral density (PSD) visualizations.
%
% Associated manuscript:
% Bertelsen N. et al.,
% "Electrophysiologically-defined excitation–inhibition autism neurosubtypes",
% Nature Neuroscience (under review).
%
% License: MIT License
% Contact: gabriele.mancini@oiit.it
% =========================================================================

clear; close all;

LFP_true = 0;  % 1 for LFP, 0 for EEG

%% -----------------------------
% Load Simulation Data
% -----------------------------

data_path = '../Data_to_plot';
file_list = dir(fullfile(data_path, '*.mat')); % or '*.csv' if CSV
if LFP_true
    load(fullfile(data_path, 'analysis_output'));
else
    load(fullfile(data_path, 'analysis_output_EEG.mat'));
end
fprintf('Data loaded');
    
%% -----------------------------
% Save figures path
% -----------------------------

if LFP_true
    save_dir = '/Users/gabrielemancini/Desktop/Bertelsen Github/In-silico/Figures/Figure 1';
else 
    save_dir = '/Users/gabrielemancini/Desktop/Bertelsen Github/In-silico/Figures/Figure 1 EEG';
end
% Create folder if it does not exist
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    fprintf('Created directory: %s\n', save_dir);
end

addpath('../../functions');

%% Parameters

labels = 0.8 * FR_exc + 0.2 * FR_inh;
thresh = 10;

%% Set default font and axis properties

set(groot, 'DefaultAxesFontSize', 10);          % applies to tick labels + axis labels
set(groot, 'DefaultTextFontSize', 10);          % applies to titles, annotations
set(groot, 'DefaultLegendFontSize', 10);        % legend text
set(groot, 'DefaultAxesFontName', 'Arial');    % axis + tick font
set(groot, 'DefaultTextFontName', 'Arial');    % text font
set(groot, 'DefaultLegendFontName', 'Arial');  % legend font
set(groot, 'DefaultAxesLineWidth', 1.5);       % axis line thickness
set(groot, 'DefaultAxesTickLabelInterpreter','tex');  % plain text ticks (no TeX)

% Extract and log-transform features
H_filter  = H_base(filter);
gamma_filter = gamma_base(filter)./max(gamma_base(filter));

nu0 = nu_0(filter);
g0 = g(filter);

% Normalize dot sizes for visibility
%dotSizes = 10 + 30 * (g0 - min(g0)) / (max(g0) - min(g0));  % sizes between 20 and 100
dotSizes = 5 + 0 * (g0 - min(g0)) / (max(g0) - min(g0));  % sizes between 20 and 100
dotSizes_d = 15 + 0 * (g0 - min(g0)) / (max(g0) - min(g0));  % sizes between 20 and 100
alphas = 0.1 + 0.8 * ( (g0 - min(g0)) / (max(g0) - min(g0)) ); 
colormap1 = (nu0 - min(nu0)) / (max(nu0) - min(nu0));
colormapnu = (nu0 - min(nu0)) / (max(nu0) - min(nu0));
colormapg = (g0 - min(g0)) / (max(g0) - min(g0));

%% K-fold 

figure('Color', 'white', 'Units', 'centimeters', 'Position', [0, 0, 8, 4]);% Labels and filter

% Labels and filter
labels = 0.8 * FR_exc + 0.2 * FR_inh;
filter = (labels < 10) & (nu_0 < 4);
labels = labels(filter);
% Extract and log-transform features
slope_filter  = slope_base(filter);
totpow_filter = Pow_base(filter);

A = { [1], [2] };

features_all = {slope_filter', totpow_filter'};
feature_names = {'Slope', 'tot pow'};
num_samples = length(labels);

k = 5; % number of folds
ci_level = 0.95;

coef_means_all = {}; % store to use later in bar plot
coef_stds_all = {};


for var = 1:2
    
    features = features_all{var};
    % Z-score normalization
    features = (features - mean(features,1)) ./ std(features,0,1);

    % Create 10-fold partition
    cv = cvpartition(num_samples, 'KFold', k);

    % Preallocate
    r_vals = zeros(k,1);
    all_predicted = nan(num_samples,1);
    all_coefs = nan(k, size(features,2)); 

    y_true_plot = [];
    y_pred_plot = [];
    regime_plot = [];
    alpha_plot = [];
    colormapnu_plot = [];

    for i = 1:k
        train_idx = training(cv, i);
        test_idx = test(cv, i);

        mdl = fitlm(features(train_idx), labels(train_idx));
        y_pred = predict(mdl, features(test_idx,:));
        y_true = labels(test_idx)';
        all_coefs(i,:) = mdl.Coefficients.Estimate(2:end)'; % exclude intercept

        % ensure they are column vectors of same size
        y_pred = y_pred(:);
        y_true = y_true(:);
        
        y_true_plot = [y_true_plot' y_true']';
        y_pred_plot = [y_pred_plot' y_pred']';
        regime_plot = [regime_plot regime(test_idx)];
        alpha_plot = [alpha_plot alphas(test_idx)]; 
        colormapnu_plot = [colormapnu_plot colormapnu(test_idx)];

        % Pearson correlation instead of R²
        r_vals(i) = corr(y_true, y_pred, 'Type', 'Pearson');


    end

    % Mean and CI of r
    r_mean = mean(r_vals);
    r_ci = prctile(r_vals, [(1-ci_level)/2*100, (1+ci_level)/2*100]);

    % Store coefficients for final bar plot
    coef_means_all{var} = mean(all_coefs,1);
    coef_stds_all{var} = std(all_coefs,0,1);

    % Scatter plot: predicted vs actual (all samples from CV)
    subplot(1, 2, A{var});
    predicted = all_predicted;
    actual = labels(:);

    hold on
    for i = 1:length(y_pred_plot)
        if regime_plot(i)==0
            scatter(y_true_plot(i), y_pred_plot(i), 3, [colormapnu_plot(i),0,0] , ...
                'filled', 'MarkerFaceAlpha', alpha_plot(i),'Marker','o');
        else
            scatter(y_true_plot(i), y_pred_plot(i), 10, [colormapnu_plot(i),0,0] , ...
                'filled', 'MarkerFaceAlpha', alpha_plot(i),'Marker','diamond');
        end
    end

    if var==1
        ylabel('Predicted FR (Hz)');
    end
    xlabel('Actual FR (Hz)'); set(gca,'FontSize',10)

    text(0.4, 0.2, { ...
        ['r = ' num2str(r_mean, '%.2f')], ...
        ['[' num2str(r_ci(1), '%.2f') ', ' num2str(r_ci(2), '%.2f') ']']}, ...
        'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');

    if var ==1
        T.R2_SLOPE(1)=r_mean;
        T.ci_SLOPE=r_ci;
    else
        T.R2_TOT(1)=r_mean;
        T.ci_TOT=r_ci;
    end
            
    title([feature_names{var}])
    h = refline(1,0);
    h.LineWidth = 1.5;
    h.LineStyle = '--';

end

save_path = fullfile(save_dir,'/panel_D.svg');
saveas(gcf, [save_path]); 

%%

% Figure setup
figure('Color','w','Units','centimeters','Position',[0 0 8 4]);

styles = {'-','--',':'};   % different line styles for g
colors = {[0.8667, 0.0863, 0.2196],[0.0863, 0.5176, 0.6235],[0. 0. 0.]}; % AMPA,GABA,LFP

g_all = -parameterdb(:,10)./parameterdb(:,13);
nu_all = parameterdb(:,53);
g_unique = unique(g_all);
nu_unique = unique(nu_all);
id = str2double(ID);
fmin = 20; 
fmax = 100; 

for i = 1:min(2, numel(nu_unique))

    subplot(1,2,i)
    hold on
    
    for j = 1:min(2, numel(g_unique))

        idx = (g_all == g_unique(20*(j-1)+1)) & (nu_all == nu_unique(i));
        k_app = find(idx, 1);
        k = find(id == k_app);
        if isempty(k), continue; end

        traces = LFP{k};
        LFP_trace  = squeeze(traces(1, trial, :));
        AMPA_trace = squeeze(traces(2, trial, :));
        GABA_trace = squeeze(traces(3, trial, :));

        [PSD_AMPA, f] = mtspectrumc(AMPA_trace, params);
        PSD_GABA      = mtspectrumc(GABA_trace, params);
        PSD_LFP       = mtspectrumc(LFP_trace, params);
 
        freq = 500*f;

        plot(freq, movmean(PSD_AMPA,200), 'Color', colors{1}, 'LineStyle', styles{j},'LineWidth',1.5)
        plot(freq, movmean(PSD_GABA,200), 'Color', colors{2}, 'LineStyle', styles{j},'LineWidth',1.5)
        plot(freq, movmean(PSD_LFP,200),  'Color', 'k', 'LineStyle', styles{j},'LineWidth',1.5)

    end
    
    set(gca,'YScale','log','XScale','log','FontSize',10,'LineWidth',1.2);
    xlim([5 150]);
    if i ==1
    ylabel('Power (a.u.)');
    end
    xlabel('Frequency (Hz)');
end

% Titles and legends
subplot(1,2,1); 
title(sprintf('\\nu = %.2f (spks/s/cell)',nu0_vals(1)));
subplot(1,2,2); 
title(sprintf('\\nu = %.2f (spks/s/cell)',nu0_vals(5)));

legendEntries = {};
for j = 1:min(3, numel(g_unique))
    legendEntries{end+1} = sprintf('AMPA g=%.2f', g_unique(7*j));
    legendEntries{end+1} = sprintf('GABA g=%.2f', g_unique(7*j));
    legendEntries{end+1} = sprintf('LFP g=%.2f',  g_unique(7*j));
end

save_path = fullfile(save_dir,'/panel_C.svg');
saveas(gcf, [save_path]); 

%% PLOT PSD 

% Unique values
nu0_vals = unique(nu0);
nColors  = 5;    % color resolution
nTraces  = 5;     % number of PSD curves per subplot
lw       = 1.5;

% Common labels
labels = 0.8*FR_exc + 0.2*FR_inh;

for i =1:length(labels)
    if (corr_base(i) < 0.009) & (isi_base(i) > 1.) & (FR_exc(i) < 2.5)
        regime(i) = 0;
    else 
        regime(i) = 1;
    end
end

% Color map: black to red
customMap = [linspace(0,1,nColors)', zeros(nColors,1), zeros(nColors,1)];

% --- First condition ---
mask1 = (labels < 10) & (nu_0 == nu0_vals(1));
PSD_1 = PSDS(mask1,:);
g_1   = g(mask1);
H_1 = H_base(mask1);
regimes_1 = regime(mask1);
% Normalize g values to [0,1] and map to colors
g1_norm = (g_1 - min(g_1)) / (max(g_1) - min(g_1));
color_idx1 = round(1 + g1_norm * (nColors - 1));

% Pick evenly spaced indices across g_1
idx1 = round(linspace(1,length(g_1),nTraces));
g_plot1 = g_1(idx1);

% --- Second condition ---
mask2 = (labels < 10) & (nu_0 == nu0_vals(5));
PSD_2 = PSDS(mask2,:);
g_2   = g(mask2);
H_2 = H_base(mask2);
regimes_2 = regime(mask2);
g2_norm = (g_2 - min(g_2)) / (max(g_2) - min(g_2));
color_idx2 = round(1 + g2_norm * (nColors - 1));

idx2 = round(linspace(1,length(g_2),nTraces));
g_plot2 = g_2(idx2);


% --- Plotting ---
% Subplot 3: true colormap (g × nu0 grid, colored by FR)
figure('Color', 'white', 'Units', 'centimeters', 'Position', [0, 0, 8, 4]);

for i = 1:nTraces
    subplot(1,2,1); hold on;
    plot(freq, movmean(PSD_1(idx1(i),:),200), ...
         'Color', customMap(color_idx1(idx1(i)),:), 'LineWidth', lw);

    subplot(1,2,2); hold on;
    plot(freq, movmean(PSD_2(idx2(i),:),200), ...
         'Color', customMap(color_idx2(idx2(i)),:), 'LineWidth', lw);
end

% Axes formatting
for k = 1:2
    subplot(1,2,k);
    set(gca, 'YScale', 'log', 'XScale', 'log', ...
             'FontSize', 10, 'LineWidth', 1.2);

    % Set log limits
    xlim([5 150]);

    % Force at least two visible log ticks
    xticks([10 100]);   % You can also try [10 30 100] if you want more
    xticklabels({'10','100'}); 

    xlabel('Frequency (Hz)')
    if k == 1
        ylabel('Power (a.u.)');
    end
end

% Titles and legends
subplot(1,2,1); 
title(sprintf('\\nu = %.2f (spks/s/cell)',nu0_vals(1)));
subplot(1,2,2); 
title(sprintf('\\nu = %.2f (spks/s/cell)',nu0_vals(5)));

save_path = fullfile(save_dir,'/panel_A.svg');
saveas(gcf, [save_path]); 

%%

figure('Color', 'white', 'Units', 'centimeters', 'Position', [0, 0, 8, 4]);

% % FOOOF settings
settings.peak_width_limits = [1, 70];  % Fixed typo
settings.max_n_peaks = 2;
settings.aperiodic_mode = 'fixed';
settings.peak_threshold =0.5;
settings.verbose = false;

f_range = f_range_slope;

for i= 1:length(g_1)
psd =PSD_1(i,:)';
fooof_results = fooof(500*f, psd, f_range, settings, false);  % 'false' to skip auto-plot
slope_1(i)=fooof_results.aperiodic_params(2);
psd =PSD_2(i,:)';
fooof_results = fooof(500*f, psd, f_range, settings, false);  % 'false' to skip auto-plot
slope_2(i)=fooof_results.aperiodic_params(2);
end

% Subplot for Group 1
subplot(2,2,1)
hold on
for i = 1:length(g_1)
    if regimes_1(i)==0
    scatter(g_1(i), -slope_1(i),5, 'k', 'filled','Marker','o');  % colormap_map is a function you'd define
    else
    scatter(g_1(i), -slope_1(i), 10, 'k','filled', 'Marker','diamond');  % colormap_map is a function you'd define
    end
end
ylabel('slope')
p1 = polyfit(g_1, -slope_1, 1); % linear fit
yfit1 = polyval(p1, g_1);
plot(g_1, yfit1, 'r-', 'LineWidth', 2)
hold off
set(gca,'FontSize',10,'LineWidth',1.2);

% Subplot for Group 2
subplot(2,2,2)
hold on
for i = 1:length(g_1)
    if regimes_2(i)==0
    scatter(g_2(i), -slope_2(i),5, 'k', 'filled','Marker','o');  % colormap_map is a function you'd define
    else
    scatter(g_2(i), -slope_2(i), 10, 'k','filled', 'Marker','diamond');  % colormap_map is a function you'd define
    end
end

ylabel('slope')
p2 = polyfit(g_2, -slope_2, 1); % linear fit
yfit2 = polyval(p2, g_2);
plot(g_2, yfit2, 'r-', 'LineWidth', 2)
hold off
set(gca, 'FontSize',10,'LineWidth',1.2);

subplot(2,2,3)
hold on
for i = 1:length(g_1)
    if regimes_1(i)==0
    scatter(g_1(i), H_1(i),5, 'k', 'filled','Marker','o');  % colormap_map is a function you'd define
    else
    scatter(g_1(i), H_1(i), 10, 'k','filled', 'Marker','diamond');  % colormap_map is a function you'd define
    end
end
xlabel('g')
ylabel('H')
p1 = polyfit(g_1, H_1, 1); % linear fit
yfit1 = polyval(p1, g_1);
plot(g_1, yfit1, 'r-', 'LineWidth', 2)
hold off
set(gca,'FontSize',10,'LineWidth',1.2);

% Subplot for Group 2
subplot(2,2,4)
hold on
for i = 1:length(g_1)
    if regimes_2(i)==0
    scatter(g_1(i), H_2(i),5, 'k', 'filled','Marker','o');  % colormap_map is a function you'd define
    else
    scatter(g_1(i), H_2(i), 10, 'k','filled', 'Marker','diamond');  % colormap_map is a function you'd define
    end
end
xlabel('g')
ylabel('H')
p2 = polyfit(g_2, H_2, 1); % linear fit
yfit2 = polyval(p2, g_2);
plot(g_2, yfit2, 'r-', 'LineWidth', 2)
hold off
set(gca,'FontSize',10,'LineWidth',1.2);

save_path = fullfile(save_dir,'/panel_B.svg');
saveas(gcf, [save_path]); 
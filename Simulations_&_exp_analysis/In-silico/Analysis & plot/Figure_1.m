% =========================================================================
%
% Author: Gabriele Mancini
% Date: 2025-12-10
%
% Description:
% Loads periodic and aperiodic spectral features extracted from simulated
% LFP/EEG data (Feature_extraction.m pipeline) together with simulation
% parameters. Performs correlation and regression analyses to relate spectral
% features to network parameters and firing rate, reproducing Fig. 1 panels.
%
% Associated manuscript:
% Bertelsen N. et al., "Electrophysiologically-defined excitation-inhibition
% autism neurosubtypes", Nature Neuroscience (under review).
%
% Figures reproduced:
% Figure 1, panels C–K
%
% License: MIT License
% Contact: gabriele.mancini@oiit.it
% =========================================================================

clear; close all;

%% -----------------------------
% Configuration
% -----------------------------
LFP_true = 0;    % 1: LFP (AMPA+GABA currents), 0: EEG proxy (Martinez et al., 2021)

%% -----------------------------
% Load data
% -----------------------------
data_path = '../Data_to_plot';
if LFP_true
    load(fullfile(data_path,'analysis_output'));
else
    load(fullfile(data_path,'analysis_output_EEG.mat'));
end
fprintf('Data loaded\n');

%% -----------------------------
% Output directory
% -----------------------------
if LFP_true
    save_dir = '/Users/gabrielemancini/Desktop/Bertelsen Github/In-silico/Figures/Figure_1_LFP';
else
    save_dir = '/Users/gabrielemancini/Desktop/Bertelsen Github/In-silico/Figures/Figure_1_EEG';
end
if ~exist(save_dir,'dir'), mkdir(save_dir); end

%% -----------------------------
% Paths and global parameters
% -----------------------------
addpath('../../functions');

labels = 0.8*FR_exc + 0.2*FR_inh;   % population firing rate (E/I weighted)
thresh = 10;                        % FR threshold (Hz)
tresh_nu = 4;                       % nu threshold (spikes/s/cell)f4

% Default figure aesthetics
set(groot,'DefaultAxesFontSize',10,'DefaultTextFontSize',10,...
    'DefaultLegendFontSize',10,'DefaultAxesFontName','Arial',...
    'DefaultTextFontName','Arial','DefaultLegendFontName','Arial',...
    'DefaultAxesLineWidth',1.5,'DefaultAxesTickLabelInterpreter','tex');

%% ========================================================================
% Panel D — FR vs g / ν slices
% ========================================================================
figure('Color','white','Units','centimeters','Position',[0 0 13 4]);

conditions = {
 @(FR_exc,FR_inh,nu_0,g) (0.8*FR_exc+0.2*FR_inh<10)&(nu_0<1.7),        1,'g','Firing Rate (Hz)','g',true;
 @(FR_exc,FR_inh,nu_0,g) (0.8*FR_exc+0.2*FR_inh<10)&(nu_0>3)&(nu_0<3.2),2,'g','',               'g',true;
 @(FR_exc,FR_inh,nu_0,g) (0.8*FR_exc+0.2*FR_inh<10)&(g>0.1)&(g<0.11),  3,'\nu (spks/s/cell)','', '\nu',false;
 @(FR_exc,FR_inh,nu_0,g) (0.8*FR_exc+0.2*FR_inh<10)&(g>0.16)&(g<0.17), 4,'\nu (spks/s/cell)','', '\nu',false;
};

titles = {'\nu=1.5','\nu=3.17','g=0.1','g=0.16'};
xlim_FR = {[0.07 0.2],[0.07 0.2],[1 5],[1 4]};
ylim_FR = {[0.2 3],[2 25],[-2 24],[-2 24]};

for i = 1:size(conditions,1)
    f = conditions{i,1}(FR_exc,FR_inh,nu_0,g);
    frE = FR_exc(f); frI = FR_inh(f);
    x = conditions{i,6}*g(f) + ~conditions{i,6}*nu_0(f);

    subplot(1,4,conditions{i,2}); hold on
    scatter(x,frE,10,'filled','MarkerFaceColor',[254 141 141]/256);
    scatter(x,frI,10,'filled','MarkerFaceColor',[113 156 248]/256);
    scatter(x,0.8*frE+0.2*frI,10,'k','filled');
    xlabel(conditions{i,3}); ylabel(conditions{i,4});
    title(titles{i}); xlim(xlim_FR{i}); ylim(ylim_FR{i});
end
saveas(gcf,fullfile(save_dir,'panel_D.png'));

%% ========================================================================
% Panel G — FR vs H and γ (regime-coded)
% ========================================================================
figure('Color','white','Units','centimeters','Position',[0 0 20 3.6]);

f = (labels<thresh)&(nu_0<tresh_nu);
labels = labels(f); nu0 = nu_0(f); g0 = g(f);
Hf = H_base(f); gamf = gamma_base(f)./max(gamma_base(f));
slopef = slope_base(f); tot_powf = Pow_base(f);
isi = isi_base(f); corrb = corr_base(f); FR_E = FR_exc(f);

% Dynamical regime classification
regime = ~(corrb<0.009 & isi>1 & FR_E<2.5);

alphas = 0.1+0.8*(g0-min(g0))/(max(g0)-min(g0));
cmap = (nu0-min(nu0))/(max(nu0)-min(nu0));

subplot(1,2,1); hold on
for i=1:numel(labels)
    scatter(labels(i),Hf(i),5+10*regime(i),[cmap(i) 0 0],'filled',...
        'Marker',char('o'*(~regime(i))+'d'*regime(i)),'MarkerFaceAlpha',alphas(i));
end
xlabel('Firing Rate (Hz)'); ylabel('H');

subplot(1,2,2); hold on
for i=1:numel(labels)
    scatter(labels(i),gamf(i),5+10*regime(i),[cmap(i) 0 0],'filled',...
        'Marker',char('o'*(~regime(i))+'d'*regime(i)),'MarkerFaceAlpha',alphas(i));
end
xlabel('Firing Rate (Hz)'); ylabel('Norm. \gamma Pow. (a.u.)');
saveas(gcf,fullfile(save_dir,'panel_G.png'));

%% ========================================================================
% Panels E–F — H, γ vs g and ν correlations
% ========================================================================

figure('Color','white','Units','centimeters','Position',[0 0 20 4]);
getR = @(x,y) corr(x(:),y(:),'Type','Pearson','Rows','complete');

% Subplot 1: H vs g
subplot(1,4,1); hold on
for i=1:numel(Hf)
    scatter(Hf(i), g0(i), 8, [cmap(i),0,0], 'filled', 'MarkerFaceAlpha', alphas(i));
end
fit = polyfit(Hf, g0, 1);
plot(sort(Hf), polyval(fit, sort(Hf)), 'k--', 'LineWidth', 1);
[r,p]=getR(Hf,g0); T.g_H=[r p];
xlabel('H'); ylabel('g'); 
text(max(Hf)*0.6, max(g0)*0.4, sprintf('r=%.2f', r), 'FontSize', 8);

% Subplot 2: H vs ν
subplot(1,4,2); hold on
for i=1:numel(Hf)
    scatter(Hf(i), nu0(i), 8, [0,0,0], 'filled', 'MarkerFaceAlpha', alphas(i));
end
fit = polyfit(Hf, nu0, 1);
plot(sort(Hf), polyval(fit, sort(Hf)), 'k--', 'LineWidth', 1);
[r,p]=getR(Hf,nu0); T.nu_H=[r p];
xlabel('H'); ylabel('\nu');
text(max(Hf)*0.6, max(nu0)*0.4, sprintf('r=%.2f', r), 'FontSize', 8);

% Subplot 3: γ vs g
subplot(1,4,3); hold on
for i=1:numel(gamf)
    scatter(gamf(i), g0(i), 8, [cmap(i),0,0], 'filled', 'MarkerFaceAlpha', alphas(i));
end
fit = polyfit(gamf, g0, 1);
plot(sort(gamf), polyval(fit, sort(gamf)), 'k--', 'LineWidth', 1);
[r,p]=getR(gamf,g0); T.g_gamma=[r p];
xlabel('\gamma'); ylabel('g');
text(max(gamf)*0.6, max(g0)*0.4, sprintf('r=%.2f', r), 'FontSize', 8);

% Subplot 4: γ vs ν
subplot(1,4,4); hold on
for i=1:numel(gamf)
    scatter(gamf(i), nu0(i), 8, [0,0,0], 'filled', 'MarkerFaceAlpha', alphas(i));
end
fit = polyfit(gamf, nu0, 1);
plot(sort(gamf), polyval(fit, sort(gamf)), 'k--', 'LineWidth', 1);
[r,p]=getR(gamf,nu0); T.nu_gamma=[r p];
xlabel('\gamma'); ylabel('\nu');
text(max(gamf)*0.6, max(nu0)*0.4, sprintf('r=%.2f', r), 'FontSize', 8);

saveas(gcf, fullfile(save_dir,'panel_E-F.png'));


%% ========================================================================
% Panel C — Dynamical regimes in (g,ν) space
% ========================================================================
load('/Users/gabrielemancini/Desktop/Bertelsen_Mancini/Simulations SFARI/colormap_regimes.mat')

figure('Color','white','Units','centimeters','Position',[0 0 6 4]);
[xq,yq]=meshgrid(linspace(min(g0),max(g0),500),linspace(min(nu0),max(nu0),500));
Z=griddata(g0,nu0,double(regime),xq,yq,'natural');
imagesc(xq(1,:),yq(:,1),Z>0.5,'AlphaData',~isnan(Z));
set(gca,'YDir','normal'); colormap(regimes); hold on
contour(xq,yq,Z,[0.5 0.5],'k','LineWidth',3)
xlabel('g'); ylabel('\nu'); title('Dynamical Regimes'); box off
saveas(gcf,fullfile(save_dir,'panel_C.png'));

%% ========================================================================
% Panels I–K — Cross-validated FR prediction
% ========================================================================

figure('Color', 'white', 'Units', 'centimeters', 'Position', [0, 0, 20, 4]);

A = { [1], [2], [3]};

features_all = {Hf', gamf', [Hf', gamf']};
feature_names = {'H', '\gamma', 'H + \gamma'};
num_samples = length(labels);

k = 5; % number of folds
ci_level = 0.95;

coef_means_all = {}; % store to use later in bar plot
coef_stds_all = {};


for var = 1:3
    
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
    
        mdl = fitlm(features(train_idx,:), labels(train_idx));
        y_pred = predict(mdl, features(test_idx,:));
        y_true = labels(test_idx)';
        all_coefs(i,:) = mdl.Coefficients.Estimate(2:end)'; % exclude intercept
    
        y_pred = y_pred(:);
        y_true = y_true(:);
    
        y_true_plot = [y_true_plot' y_true']';
        y_pred_plot = [y_pred_plot' y_pred']';
        regime_plot = [regime_plot regime(test_idx)];
        alpha_plot = [alpha_plot alphas(test_idx)]; 
        colormapnu_plot = [colormapnu_plot cmap(test_idx)];
    
        % Pearson correlation instead of R²
        r_vals(i) = corr(y_true, y_pred, 'Type', 'Pearson');

    end
    
    % Mean and CI of r
    r_mean = mean(r_vals);
    r_ci = prctile(r_vals, [(1-ci_level)/2*100, (1+ci_level)/2*100]);

    % Store coefficients for final bar plot
    coef_means_all{var} = mean(all_coefs,1);
    coef_stds_all{var} = std(all_coefs,0,1);

    if var == 3
    B_c = mean(all_coefs,1);
    save('/Users/gabrielemancini/Desktop/Bertelsen_Mancini/experimental_data/DATA/code/Weights_simulations.mat', 'B_c');
    end

    % Scatter plot: predicted vs actual (all samples from CV)
    subplot(1, 4, A{var});
    predicted = all_predicted;
    actual = labels(:);

    hold on
    for i = 1:length(y_pred_plot)
        if regime_plot(i)==0
            scatter(y_true_plot(i), y_pred_plot(i), 5, [colormapnu_plot(i),0,0] , ...
                'filled', 'MarkerFaceAlpha', alpha_plot(i),'Marker','o');
        else
            scatter(y_true_plot(i), y_pred_plot(i), 15, [colormapnu_plot(i),0,0] , ...
                'filled', 'MarkerFaceAlpha', alpha_plot(i),'Marker','diamond');
        end
    end

    if var==1
        ylabel('Predicted FR (Hz)');
    end
    xlabel('Actual FR (Hz)'); set(gca,'FontSize',10)

    % Annotate with R² mean and CI
    text(0.2, 0.1, ...
        ['r = ' num2str(r_mean, '%.2f') ...
        ' [' num2str(r_ci(1), '%.2f') ', ' num2str(r_ci(2), '%.2f') ']'], ...
        'Units', 'normalized', 'FontSize', 8, 'FontWeight', 'bold')

    if var ==1
        T.R2_H(1) = r_mean;
        T.ci_H = r_ci;
    elseif var ==2
        T.R2_gamma(1) = r_mean;
        T.ci_gamma = r_ci;
    elseif var ==2 
        T.R2_HGAMMA(1) = r_mean;
        T.ci_HGAMMA = r_ci;
    elseif var ==3 
        T.R2_slope(1) = r_mean;
        T.ci_slope = r_ci;
    elseif var ==4 
        T.R2_power(1) = r_mean;
        T.ci_power = r_ci;
    end

    title([feature_names{var}])
    h = refline(1,0);
    h.LineWidth = 1.5;
    h.LineStyle = '--';

end

% === Bar plot of mean weights from the combined (H+gamma) model ===
subplot(1, 4, 4);
coef_mean = coef_means_all{3};
coef_std = coef_stds_all{3};

bar(coef_mean, 'FaceColor', [0.5, 0.5, 0.5]); hold on
for j = 1:length(coef_mean)
    plot([j j], [coef_mean(j)-coef_std(j), coef_mean(j)+coef_std(j)], ...
        'k', 'LineWidth', 1.5);
end
xticks(1:length(coef_mean))
xticklabels({'H','\gamma'});
box('off')
ylabel('Coefficient'); set(gca,'FontSize',8)

save_path = fullfile(save_dir,'/panel_H-I.png');
saveas(gcf, [save_path]); 


%% ========================================================================
% Compute performance for other models
% ========================================================================

features_all = {slopef', tot_powf'};
feature_names = {'Slope','tot pow'};
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

    for i = 1:k
        train_idx = training(cv, i);
        test_idx = test(cv, i);
    
        mdl = fitlm(features(train_idx,:), labels(train_idx));
        y_pred = predict(mdl, features(test_idx,:));
        y_true = labels(test_idx)';
        all_coefs(i,:) = mdl.Coefficients.Estimate(2:end)'; % exclude intercept
    
        y_pred = y_pred(:);
        y_true = y_true(:);
   
        % Pearson correlation instead of R²
        r_vals(i) = corr(y_true, y_pred, 'Type', 'Pearson');

    end
    
    % Mean and CI of r
    r_mean = mean(r_vals);
    r_ci = prctile(r_vals, [(1-ci_level)/2*100, (1+ci_level)/2*100]);


   if var ==1 
        T.R2_slope(1) = r_mean;
        T.ci_slope = r_ci;
   else
        T.R2_power(1) = r_mean;
        T.ci_power = r_ci;
    end


end

%% -----------------------------
% Display summary statistics
% -----------------------------
disp(T)

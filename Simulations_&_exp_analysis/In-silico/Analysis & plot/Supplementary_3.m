% =========================================================================
%  SUPPLEMENTARY FIGURE 2 — RESTING POTENTIAL MANIPULATION
%
% Author: Gabriele Mancini
% Date: 2025-12-10
%
% Description:
%  Loads simulation data, extracts features, averages trials, and plots
%  Δ Firing rate, Δ Hurst, and Δ γ power as a function of resting potential.
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


%% Loop over the manipulations

for camk=0:2

% ------------------------------------------------------------------------
%  USER SETTINGS
% -------------------------------------------------------------------------
trials      = 2;

base_path   = '/Users/gabrielemancini/Desktop/Bertelsen Github/In-silico/Simulations/Simulation results';
save_dir    = '/Users/gabrielemancini/Desktop/Bertelsen Github/In-silico/Figures/Supplementaries/Supplementary 3';

% Create folder if it does not exist
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    fprintf('Created directory: %s\n', save_dir);
end

addpath('../../functions');

%% ------------------------------------------------------------------------
%  SELECT EXPERIMENT PATH
% -------------------------------------------------------------------------
switch camk
    case 1, root_path = fullfile(base_path, 'camk');
    case 0, root_path = fullfile(base_path, 'hm4di');
    otherwise, root_path = fullfile(base_path, 'PV');
end

path_sim    = fullfile(root_path, 'data');
path_stats  = fullfile(root_path, 'metadata/stats_db.csv');
path_params = fullfile(root_path, 'metadata/parameter_db.csv');

addpath('../functions');

%% ------------------------------------------------------------------------
%  LOAD SIMULATION DATA
% -------------------------------------------------------------------------
files = dir(fullfile(path_sim, '*.npy'));
files = {files.name};

LFP = cell(1, numel(files));
ID  = cell(1, numel(files));

for k = 1:numel(files)
    fprintf('\rLoading LFP file %d / %d (%.0f%%)', k, numel(files), 100 * k / numel(files));
    LFP{k} = readNPY(fullfile(path_sim, files{k}));
    ID{k}  = extractAfter(erase(files{k}, '.npy'), 'ID');
end
fprintf('\n');

%% Load metadata
statsdb     = readmatrix(path_stats);   statsdb(1,:) = [];
parameterdb = readmatrix(path_params);  parameterdb(1,:) = []; parameterdb(:,1:2) = [];

%% ------------------------------------------------------------------------
%  PARAMETERS
% -------------------------------------------------------------------------
fs               = 500;
params.tapers    = [3 5];
f_range          = [1 80];
f_range_gamma    = [30 100];
f_range_slope    = [1 150];
Hparams          = get_hurst_parameters;

%% ------------------------------------------------------------------------
%  SORT DATA BY ID
% -------------------------------------------------------------------------
ID_num = str2double(ID);
[~, idx] = sort(ID_num);

ID_sorted  = ID(idx);
LFP_sorted = LFP(idx);

%% ------------------------------------------------------------------------
%  FEATURE EXTRACTION
% -------------------------------------------------------------------------
count = 1;

for i = 1:numel(LFP)

    id      = str2double(ID_sorted{i});
    nu_val  = parameterdb(parameterdb(:,59) == id, 53);

    % Skip pathological ν values
    if nu_val >= 10
        continue;
    end

    for tr = 1:trials

        LFP_trace = squeeze(LFP_sorted{i}(1, tr, :));
        LFP_detr  = detrend(LFP_trace);

        % PSD & Power
        [tot_pow, ~, ~] = compute_power(LFP_trace, fs, f_range);
        tot_power(tr,count) = tot_pow;

        % Firing rate
        fr_exc_app = statsdb(statsdb(:,10)==id,1); 
        fr_inh_app = statsdb(statsdb(:,10)==id,5); 
        fr_exc(tr,count) = fr_exc_app(1); 
        fr_inh(tr,count) = fr_inh_app(1);

        % Synaptic parameters
        idx_param = (parameterdb(:,59) == id);
        g_param(tr,count)  = -parameterdb(idx_param,10)./parameterdb(idx_param,13);
        g_e(tr,count)      =  parameterdb(idx_param,10);
        g_i(tr,count)      = -parameterdb(idx_param,13);
        EL(tr,count)       =  parameterdb(idx_param,39);
        g_EE(tr,count)     =  parameterdb(idx_param,12);
        nu(tr,count)       =  nu_val;

        % Hurst
        H_LFP(tr,count) = bfn_mfin_ml(LFP_detr, 'filter', Hparams.filter, ...
                                      'lb',Hparams.lb,'ub',Hparams.ub, ...
                                      'verbose',Hparams.verbose);

        % γ Power and slope (FOOOF)
        gamma_pow_LFP(tr,count) = compute_psd_area(LFP_detr, fs, f_range_gamma);
        [a,f0,ampl,bw] = compute_spectral_slope_fooof(LFP_trace, fs, f_range_slope);
        slope_LFP_fooof(tr,count)    = a;
        gamma_freq_LFP_fooof(tr,count) = f0;
        amplitude_LFP_fooof(tr,count) = ampl;
        bw_LFP_fooof(tr,count)        = bw;

        % Multitaper PSD & slope
        [S, f] = mtspectrumc(LFP_trace, params);
        PSDs(tr,count,:) = S';
        slope(tr,count)  = compute_spectral_slope(LFP_trace, fs, f_range_slope);

        count = count + 1;
    end
end

%% ------------------------------------------------------------------------
%  AVERAGE TRIALS
% -------------------------------------------------------------------------
avg = @(x) (x(1,1:2:end-1) + x(2,2:2:end)) / 2;

FR_exc         = avg(fr_exc);
FR_inh         = avg(fr_inh);
FR_combined    = 0.8 * FR_exc + 0.2 * FR_inh;

g              = avg(g_param);
g_E            = avg(g_EE);
EL_mean        = avg(EL);
nu_0           = avg(nu);

H_base         = avg(H_LFP);
gamma_base     = avg(gamma_pow_LFP);
slope_base     = avg(slope_LFP_fooof);
Pow_base       = avg(tot_power);

gamma_freq_base     = avg(gamma_freq_LFP_fooof);
gamma_amplitude_base = avg(amplitude_LFP_fooof);
bw_base             = avg(bw_LFP_fooof);

PSDs_mean = avg(PSDs);
PSDS      = squeeze(PSDs_mean);

%% ------------------------------------------------------------------------
%  SELECT MANIPULATED PARAMETER (RESTING POTENTIAL)
% -------------------------------------------------------------------------
switch camk
    case 1, param = parameterdb(:,25);
    case 0, param = parameterdb(:,25);
    otherwise, param = parameterdb(:,39);
end

%% ------------------------------------------------------------------------
%  NORMALISE TO CONTROL (Δ values)
% -------------------------------------------------------------------------
FR_combined_pc = FR_combined - FR_combined(1);
H_base_pc      = H_base      - H_base(1);
gamma_base_pc  = gamma_base  - gamma_base(1);

%% ------------------------------------------------------------------------
%  PLOTTING SETTINGS
% -------------------------------------------------------------------------
set(groot,'DefaultAxesFontSize',10, 'DefaultTextFontSize',10,...
          'DefaultAxesFontName','Arial','DefaultTextFontName','Arial',...
          'DefaultAxesLineWidth',1.5);

figure_width_cm  = 4;
figure_height_cm = 3;

colormap_vals = (g_E - min(g_E)) ./ (max(g_E) - min(g_E));
if camk == 1
    colormap_vals = zeros(size(colormap_vals));
end

%% ========================================================================
%  PANEL A — Δ Firing Rate
% ========================================================================
figure('Color','white','Units','centimeters',...
       'Position',[0 0 figure_width_cm figure_height_cm]); hold on;

for k = 1:numel(FR_combined_pc)-5
    scatter(param(k), FR_combined_pc(k), 20, repmat(colormap_vals(k),1,3), 'filled');
end

xline(-70,'b--','LineWidth',1.5);
xlabel('Resting Potential (mV)');
ylabel('\Delta Firing Rate (Hz)');

set(gca,'XDir','reverse');

if camk ==0
saveas(gcf, fullfile(save_dir,'hm4di_FR.svg'));
elseif camk ==1
saveas(gcf, fullfile(save_dir,'camk_FR.svg'));
else
saveas(gcf, fullfile(save_dir,'pv_FR.svg'));
end


%% ========================================================================
%  PANEL B — Δ Hurst
% ========================================================================
figure('Color','white','Units','centimeters',...
       'Position',[0 0 figure_width_cm figure_height_cm]); hold on;

for k = 1:numel(H_base_pc)-5
    scatter(param(k), H_base_pc(k), 20, repmat(colormap_vals(k),1,3), 'filled');
end

xline(-70,'b--','LineWidth',1.5);
xlabel('Resting Potential (mV)');
ylabel('\Delta H');

set(gca,'XDir','reverse');

if camk ==0

saveas(gcf, fullfile(save_dir,'hm4di_H.svg'));
elseif camk ==1
saveas(gcf, fullfile(save_dir,'camk_H.svg'));
else
saveas(gcf, fullfile(save_dir,'pv_H.svg'));
end

%% ========================================================================
%  PANEL C — Δ γ Power
% ========================================================================
figure('Color','white','Units','centimeters',...
       'Position',[0 0 figure_width_cm figure_height_cm]); hold on;

for k = 1:numel(gamma_base_pc)-5
    scatter(param(k), gamma_base_pc(k), 20, repmat(colormap_vals(k),1,3), 'filled');
end

xline(-70,'b--','LineWidth',1.5);
xlabel('Resting Potential (mV)');
ylabel('\Delta \gamma Power (a.u.)');

set(gca,'XDir','reverse');
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))

if camk ==0
saveas(gcf, fullfile(save_dir,'hm4di_gamma.svg'));
elseif camk ==1
saveas(gcf, fullfile(save_dir,'camk_gamma.svg'));
else
saveas(gcf, fullfile(save_dir,'pv_gamma.svg'));
end

end

%% ========================================================================
fprintf('Finished generating Supplementary Figure 2.\n');

clear; close all;


trials = 2;

base_path = '../Simulations/Simulation results';

save_dir    = '../Figures/Supplementaries/Supplementary NI';

% Create folder if it does not exist
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    fprintf('Created directory: %s\n', save_dir);
end


root_path = fullfile(base_path, 'data_NI');

path_stats = fullfile(root_path, 'metadata/stats_db.csv');
path_params = fullfile(root_path, 'metadata/parameter_db.csv');

addpath('../../functions');

files = dir(fullfile(root_path, '*.npy'));
files = {files.name};
LFP = cell(1, numel(files));
ID = cell(1, numel(files));

for k = 1:numel(files)
    fprintf('\rLoading LFP data file %d (%.0f%%)', k, 100 * k / numel(files));
    LFP{k} = readNPY(fullfile(root_path, files{k}));
    ID{k} = extractAfter(erase(files{k}, '.npy'), 'ID');
end


statsdb = readmatrix(path_stats); statsdb = statsdb(1:end, :);
parameterdb = readmatrix(path_params); parameterdb = parameterdb(1:end, 3:end);
parameterdb(1,:) = [];
statsdb(1,:) = [];



%% Parameters

fs = 500;
f_range = [1,80]; f_range_gamma = [30,100]; f_range_slope = [1,150];
params.tapers = [3 5];
Hparams = get_hurst_parameters;

%% Compute Features

count = 1;

% Convert to numeric
ID_num = str2double(ID);
% Get sorting indices based on numeric values
[~, idx] = sort(ID_num);
% Apply sort to original cell array
ID_sorted = ID(idx);
LFP_sorted =LFP(idx);

for i = 1:numel(LFP)

    id = str2double(ID_sorted{i});

    nu_val = parameterdb(parameterdb(:, 59) == id, 53);

    if nu_val >= 10, continue; end

    for trial = 1:2

        LFP_trace = squeeze(LFP_sorted{i}(1, trial, :));

        % Compute PSD power
        [tot_pow, gamma_pow, ~] = compute_power(LFP_trace, fs, f_range);
        tot_power(trial,count) = tot_pow;


        fr_exc_app = statsdb(statsdb(:,10)==id,1);
        fr_inh_app = statsdb(statsdb(:,10)==id,5);

        fr_exc(trial,count) = fr_exc_app(1);
        fr_inh(trial,count) = fr_inh_app(1);


        % Compute g and nu
        param_idx = parameterdb(:,59) == id;

        g_param(trial,count) = -parameterdb(parameterdb(:,59)==id,10)./parameterdb(parameterdb(:,59)==id,13);
        
        EL(trial,count)  = parameterdb(parameterdb(:,59)==id,39);
        g_EE(trial,count) = parameterdb(parameterdb(:,59)==id,12);

        nu(trial,count) = nu_val;

        %Hurst & spectral slope (if needed)

        [H_LFP(trial,count),~,~,~] = bfn_mfin_ml(detrend(LFP_trace), 'filter', Hparams.filter, 'lb',Hparams.lb, 'ub',Hparams.ub, 'verbose',Hparams.verbose);
        gamma_pow_LFP(trial,count) = compute_psd_area(detrend(LFP_trace), fs, f_range_gamma);
        [a,f,ampl,bw] = compute_spectral_slope_fooof(LFP_trace, fs, f_range_slope);
        slope_LFP_fooof(trial,count) = a; gamma_freq_LFP_fooof(trial,count) = f;
        amplitude_LFP_fooof(trial,count) = ampl; bw_LFP_fooof(trial,count) = bw;

        [S, f] = mtspectrumc(LFP_trace, params);
        PSDs(trial,count,:) = S';
        slope(trial,count) = compute_spectral_slope(LFP_trace, fs, f_range_slope);

        count = count + 1;

    end

end


%% Average Trials

avg_trials = @(x) (x(1,1:2:end-1) + x(2,2:2:end)) / 2;
FR_exc = avg_trials(fr_exc); FR_inh = avg_trials(fr_inh);
g = avg_trials(g_param); nu_0 = avg_trials(nu);
EL_mean = avg_trials(EL); nu_0 = avg_trials(nu);
g_E = avg_trials(g_EE);
H_base = avg_trials(H_LFP); gamma_base = avg_trials(gamma_pow_LFP);
slope_base = avg_trials(slope_LFP_fooof); Pow_base = avg_trials(tot_power);
gamma_freq_base = avg_trials(gamma_freq_LFP_fooof);
gamma_amplitude_base = avg_trials(amplitude_LFP_fooof);
bw_base = avg_trials(bw_LFP_fooof);
PSDs_mean =  (PSDs(1,1:2:end-1,:) + PSDs(2,2:2:end,:)) /2;
PSDS = squeeze(PSDs_mean);


%%

numb = size(PSDS,1);
coloramp_psd = parula(numb);
FR_combined = 0.2 * FR_inh + 0.8 * FR_exc;



%%


set(groot, 'DefaultAxesFontSize', 10);          % applies to tick labels + axis labels
set(groot, 'DefaultTextFontSize', 10);          % applies to titles, annotations
set(groot, 'DefaultLegendFontSize', 10);        % legend text
set(groot, 'DefaultAxesFontName', 'Arial');    % axis + tick font
set(groot, 'DefaultTextFontName', 'Arial');    % text font
set(groot, 'DefaultLegendFontName', 'Arial');  % legend font
set(groot, 'DefaultAxesLineWidth', 1.5);       % axis line thickness
set(groot, 'DefaultAxesTickLabelInterpreter','tex');  % plain text ticks (no TeX)


param = -parameterdb(:,10)./parameterdb(:,13);

numb = size(PSDS,1);
coloramp_psd = parula(numb);
FR_combined = 0.2 * FR_inh + 0.8 * FR_exc;

% % --- Convert to percentage change relative to first element ---
FR_exc_pc      = 100 * (FR_exc      - FR_exc(1))      / FR_exc(1);
FR_inh_pc      = 100 * (FR_inh      - FR_inh(1))      / FR_inh(1);
FR_combined_pc = 100 * (FR_combined - FR_combined(1)) / FR_combined(1);
H_base_pc      = 100 * (H_base      - H_base(1))      / H_base(1);
gamma_base_pc  = 100 * (gamma_base  - gamma_base(1))  / gamma_base(1);
           
color_AMPA = [0.8667, 0.0863, 0.2196];
color_GABA = [0.0863, 0.5176, 0.6235];


%%

NI = 1000;
scal = [0.3,0.5,1,2,3];

figure_width_cm = 16.;
figure_height_cm = 4.;

% Create figure
fig = figure;
set(fig, 'Color', 'white');
set(fig, 'Units', 'centimeters')
set(fig, 'Position', [0, 0, figure_width_cm, figure_height_cm]);


for i =1:5
subplot(1,5,i)
plot(param(4*(i-1)+1:4*i), FR_exc(4*(i-1)+1:4*i),'color',color_AMPA,'Marker','o','LineWidth',1.5,'MarkerFaceColor',color_AMPA,'MarkerSize',3)
hold on
plot(param(4*(i-1)+1:4*i), FR_inh(4*(i-1)+1:4*i), 'color',color_GABA,'Marker','o','LineWidth',1.5,'MarkerFaceColor',color_GABA,'MarkerSize',3)
plot(param(4*(i-1)+1:4*i), FR_combined(4*(i-1)+1:4*i), 'color', 'k','Marker','o','LineWidth',1.5,'MarkerFaceColor','k','MarkerSize',3)
xlabel('g')
title(['N_{inh} = ' num2str(int32(NI / scal(i)))]);
if i==1
    ylabel("Firing Rate (Hz)")
end
box('off')
end

saveas(gcf, fullfile(save_dir,'Supplementary_NI.svg'));
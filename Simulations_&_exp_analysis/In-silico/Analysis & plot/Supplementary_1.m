clear all
close all

addpath functions/ functions/chronux_2_12/

LFP_ratio =readmatrix('../Simulations/Simulation results/Independent_populations/LFP_cost_ratio.csv');
LFP_sum =readmatrix('../Simulations/Simulation results/Independent_populations/LFP_cost_sum.csv');
LFP_g =readmatrix('../Simulations/Simulation results/Independent_populations/LFP_g.csv');

LFP_all{1} = LFP_ratio(2:226,500:end);
LFP_all{2} = LFP_sum(2:226,500:end);
LFP_all{3} = LFP_g(2:226,500:end);

param{1} =readmatrix('../Simulations/Simulation results/Independent_populations/results_cost_ratio.csv');
param{2} =readmatrix('../Simulations/Simulation results/Independent_populations/results_cost_sum.csv');
param{3} =readmatrix('../Simulations/Simulation results/Independent_populations/results_g.csv');

SAVE_DIR  = '../Figures/Figure 2';

%% compute hurst



Hparams = get_hurst_parameters;

for i=1:225

LFP=LFP_all{1}(i,:);
[H, nfcor, fcor, att] = bfn_mfin_ml(detrend(LFP'), 'filter', Hparams.filter, ...
                                            'lb',Hparams.lb, ...
                                            'ub',Hparams.ub, ...
                                            'verbose',Hparams.verbose);
H_1(i)=H;


LFP=LFP_all{2}(i,:);
[H, nfcor, fcor, att] = bfn_mfin_ml(detrend(LFP'), 'filter', Hparams.filter, ...
                                            'lb',Hparams.lb, ...
                                            'ub',Hparams.ub, ...
                                            'verbose',Hparams.verbose);
H_2(i)=H;


LFP=LFP_all{3}(i,:);
[H, nfcor, fcor, att] = bfn_mfin_ml(detrend(LFP'), 'filter', Hparams.filter, ...
                                            'lb',Hparams.lb, ...
                                            'ub',Hparams.ub, ...
                                            'verbose',Hparams.verbose);
H_3(i)=H;

end


param{1}(:,end+1)=H_1';
param{2}(:,end+1)=H_2';
param{3}(:,end+1)=H_3';


%% COmpute mean across trials

trials=15;
for i=1:(size(param{1},1)/trials)
data_mu(i,:)=mean(param{1}((i-1)*trials+1:trials*i,:));
data_std(i,:)=std(param{1}((i-1)*trials+1:trials*i,:));
end

for i=1:(size(param{2},1)/trials)
data2_mu(i,:)=mean(param{2}((i-1)*trials+1:trials*i,:));
data2_std(i,:)=std(param{2}((i-1)*trials+1:trials*i,:));
end

for i=1:(size(param{3},1)/trials)
data3_mu(i,:)=mean(param{3}((i-1)*trials+1:trials*i,:));
data3_std(i,:)=std(param{3}((i-1)*trials+1:trials*i,:));
end

data3_mu(:,1)=abs(data3_mu(:,1));

%% Compute correlations


[rho_1,pval_1] = corr(param{1}(:,1), param{1}(:,2));
[rho_2,pval_2] = corr(param{2}(:,1), param{2}(:,2));
[rho_3,pval_3] = corr(-param{3}(:,1), param{3}(:,2));


[rho_1_H,pval_1_H] = corr(param{1}(:,1), param{1}(:,3));
[rho_2_H,pval_2_H] = corr(param{2}(:,1), param{2}(:,3));
[rho_3_H,pval_3_H] = corr(-param{3}(:,1), param{3}(:,3));



%% PLOT

figure_width_cm = 18.;
figure_height_cm = 12.;
font_name = 'Arial';
font_size = 7;
font_size_title = 8;

% Specify font name and size
fig = figure;
set(fig, 'Color', 'white');
set(fig, 'Units', 'centimeters');
set(fig, 'Position', [0, 0, figure_width_cm, figure_height_cm]);

subplot(2,3,1)
% Plot mean data
plot(data_mu(:,1), data_mu(:,2), 'b-', 'LineWidth', 2.);
hold on;
%xline(2.6, '--k', 'LineWidth', 2.)
fill([data_mu(:,1); flipud(data_mu(:,1))], [data_mu(:,2) + data_std(:,2); flipud(data_mu(:,2) - data_std(:,2))], 'b', 'FaceAlpha', 0.2);
xlabel('Total Firing Rate (Hz)', 'FontSize', font_size, 'FontName', font_name);
ylabel('1/f slope', 'FontSize', font_size, 'FontName', font_name);
title('FR_E / FR_I = cost', 'FontSize', font_size_title, 'FontName', font_name);
ylim([-1.5,-0.1])
text_x=xlim;
text_y=ylim;
text(text_x(2)-0.2*(text_x(2)-text_x(1)),text_y(2)-.8*(text_y(2)-text_y(1)), ['r = ' sprintf('%.2f', rho_1)], 'FontWeight', 'bold', 'Color', 'black', 'FontSize', font_size, 'FontName', font_name, 'HorizontalAlignment', 'center');
text(text_x(2)-0.2*(text_x(2)-text_x(1)),text_y(2)-.9*(text_y(2)-text_y(1)), ['p = ' sprintf('%.3f', pval_1)], 'FontWeight', 'bold', 'Color', 'black', 'FontSize', font_size, 'FontName', font_name, 'HorizontalAlignment', 'center');
hold off;
set(gca, 'FontName', font_name, 'FontSize', font_size);

subplot(2,3,2)
% Plot mean data
plot(data2_mu(:,1), data2_mu(:,2), 'b-', 'LineWidth', 2.);
hold on;
%xline(2/5, '--k', 'LineWidth', 2.)
fill([data2_mu(:,1); flipud(data2_mu(:,1))], [data2_mu(:,2) + data2_std(:,2); flipud(data2_mu(:,2) - data2_std(:,2))], 'b', 'FaceAlpha', 0.2);
xlabel('FR_E / FR_I', 'FontSize', font_size, 'FontName', font_name);
title('Total Firing Rate = cost', 'FontSize', font_size_title, 'FontName', font_name);
xlim([0,0.5])
ylim([-1.5,-0.1])
text_x=xlim;
text_y=ylim;
text(text_x(2)-0.2*(text_x(2)-text_x(1)),text_y(2)-.8*(text_y(2)-text_y(1)), ['r = ' sprintf('%.2f', rho_2)], 'FontWeight', 'bold', 'Color', 'black', 'FontSize', font_size, 'FontName', font_name, 'HorizontalAlignment', 'center');
text(text_x(2)-0.2*(text_x(2)-text_x(1)),text_y(2)-.9*(text_y(2)-text_y(1)), ['p = ' sprintf('%.3f', pval_2)], 'FontWeight', 'bold', 'Color', 'black', 'FontSize', font_size, 'FontName', font_name, 'HorizontalAlignment', 'center');

hold off;
set(gca, 'FontName', font_name, 'FontSize', font_size);

subplot(2,3,3)
% Plot mean data
plot(data3_mu(:,1), data3_mu(:,2), 'b-', 'LineWidth', 2.);
hold on;
%xline(0.2/1.5, '--k', 'LineWidth', 2.)
fill([data3_mu(:,1); flipud(data3_mu(:,1))], [data3_mu(:,2) + data3_std(:,2); flipud(data3_mu(:,2) - data3_std(:,2))], 'b', 'FaceAlpha', 0.2);
xlabel('g_E / g_I', 'FontSize', font_size, 'FontName', font_name);
title('Total Firing Rate = cost', 'FontSize', font_size_title, 'FontName', font_name);
xlim([0.05,0.52])
ylim([-1.5,-0.1])
text_x=xlim;
text_y=ylim;
text(text_x(2)-0.2*(text_x(2)-text_x(1)),text_y(2)-.8*(text_y(2)-text_y(1)), ['r = ' sprintf('%.2f', rho_3)], 'FontWeight', 'bold', 'Color', 'black', 'FontSize', font_size, 'FontName', font_name, 'HorizontalAlignment', 'center');
text(text_x(2)-0.2*(text_x(2)-text_x(1)),text_y(2)-.9*(text_y(2)-text_y(1)), ['p = ' sprintf('%.3f', pval_3)], 'FontWeight', 'bold', 'Color', 'black', 'FontSize', font_size, 'FontName', font_name, 'HorizontalAlignment', 'center');
hold off;
set(gca, 'FontName', font_name, 'FontSize', font_size);


subplot(2,3,4)
% Plot mean data
plot(data_mu(:,1), data_mu(:,3), 'b-', 'LineWidth', 2.);
hold on;
%xline(2.6, '--k', 'LineWidth', 2.)
fill([data_mu(:,1); flipud(data_mu(:,1))], [data_mu(:,3) + data_std(:,3); flipud(data_mu(:,3) - data_std(:,3))], 'b', 'FaceAlpha', 0.2);
xlabel('Total Firing Rate (Hz)', 'FontSize', font_size, 'FontName', font_name);
ylabel('H', 'FontSize', font_size, 'FontName', font_name);
title('FR_E / FR_I = cost', 'FontSize', font_size_title, 'FontName', font_name);
hold off;
xlim([2.2,8.])
ylim([0.92,1.2])
text_x=xlim;
text_y=ylim;
text(text_x(2)-0.2*(text_x(2)-text_x(1)),text_y(2)-.1*(text_y(2)-text_y(1)), ['r = ' sprintf('%.2f', rho_1_H)], 'FontWeight', 'bold', 'Color', 'black', 'FontSize', font_size, 'FontName', font_name, 'HorizontalAlignment', 'center');
text(text_x(2)-0.2*(text_x(2)-text_x(1)),text_y(2)-.2*(text_y(2)-text_y(1)), ['p = ' sprintf('%.3f', pval_1_H)], 'FontWeight', 'bold', 'Color', 'black', 'FontSize', font_size, 'FontName', font_name, 'HorizontalAlignment', 'center');
set(gca, 'FontName', font_name, 'FontSize', font_size);


subplot(2,3,5)
% Plot mean data
plot(data2_mu(:,1), data2_mu(:,3), 'b-', 'LineWidth', 2.);
hold on;
%xline(2/5, '--k', 'LineWidth', 2.)
fill([data2_mu(:,1); flipud(data2_mu(:,1))], [data2_mu(:,3) + data2_std(:,3); flipud(data2_mu(:,3) - data2_std(:,3))], 'b', 'FaceAlpha', 0.2);
xlabel('FR_E / FR_I', 'FontSize', font_size, 'FontName', font_name);
title('Total Firing Rate = cost', 'FontSize', font_size_title, 'FontName', font_name);
xlim([0.,0.5])
ylim([0.92,1.2])
text_x=xlim;
text_y=ylim;
text(text_x(2)-0.2*(text_x(2)-text_x(1)),text_y(2)-.1*(text_y(2)-text_y(1)), ['r = ' sprintf('%.2f', rho_2_H)], 'FontWeight', 'bold', 'Color', 'black', 'FontSize', font_size, 'FontName', font_name, 'HorizontalAlignment', 'center');
text(text_x(2)-0.2*(text_x(2)-text_x(1)),text_y(2)-.2*(text_y(2)-text_y(1)), ['p = ' sprintf('%.3f', pval_2_H)], 'FontWeight', 'bold', 'Color', 'black', 'FontSize', font_size, 'FontName', font_name, 'HorizontalAlignment', 'center');
hold off;
set(gca, 'FontName', font_name, 'FontSize', font_size);

subplot(2,3,6)
% Plot mean data
plot(data3_mu(:,1), data3_mu(:,3), 'b-', 'LineWidth', 2.);
hold on;
%xline(0.2/1.5, '--k', 'LineWidth', 2.)
fill([data3_mu(:,1); flipud(data3_mu(:,1))], [data3_mu(:,3) + data3_std(:,3); flipud(data3_mu(:,3) - data3_std(:,3))], 'b', 'FaceAlpha', 0.2);
xlabel('g_E / g_I', 'FontSize', font_size, 'FontName', font_name);
title('Total Firing Rate = cost', 'FontSize', font_size_title, 'FontName', font_name);
xlim([0.,0.52])
ylim([0.92,1.2])
text_x=xlim;
text_y=ylim;
text(text_x(2)-0.2*(text_x(2)-text_x(1)),text_y(2)-.1*(text_y(2)-text_y(1)), ['r = ' sprintf('%.2f', rho_3_H)], 'FontWeight', 'bold', 'Color', 'black', 'FontSize', font_size, 'FontName', font_name, 'HorizontalAlignment', 'center');
text(text_x(2)-0.2*(text_x(2)-text_x(1)),text_y(2)-.2*(text_y(2)-text_y(1)), ['p = ' sprintf('%.3f', pval_3_H)], 'FontWeight', 'bold', 'Color', 'black', 'FontSize', font_size, 'FontName', font_name, 'HorizontalAlignment', 'center');
hold off;
set(gca, 'FontName', font_name, 'FontSize', font_size);




saveas(fig_d, fullfile(SAVE_DIR, 'Supplementary_1.svg'));

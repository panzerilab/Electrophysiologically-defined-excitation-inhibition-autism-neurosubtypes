function slope_stats_out = perform_slope_stats(results_all,slope)
% Function to do slope stats
% Inputs: 
% Results_all: a cell array of length of n_subjects iside each cell there 
% is a structure with the following fields:
% - type: the type of the subject (exp, ctrl,...)
% - params: a structure containing the parameters of all the analyses. for
% example params.slope contains the parameters of the slope analyses
% - results: a structure containing the results of various analyses. for
% example results.slope contains the results of the slope analyses
% slope: a structure containing the parameters for the group level 
% statistical analyses with the following fields:
% -frequency: an array lenght two of positive numbers ([f1,f2], f2>f1) which determines the beginiing and ending of the range that you have fitted the line
% Outputs:
% slope_stats_out: a structure containing various fields each with the
% results of the statistical analyses performed:
% - spectrum_stats_out.(phase_name).(region_name).(type_name), is a structure
% containing various info for the data in phase = phase_name
% (baseline,cno,..) and the region_name (PFC,rs,...) for the data type
% (exp,ctrl,...) with the following fields:
%   - spec_freq: the frequency at which the spectrum was calculated
%   - spec_mean: the average spectrum in log10
%   - spec_se: the standard error of the average of spectrum in log10 scale
%   - beta: the average of intercept and slope values for data of that type
% - spectrum_stats_out.(phase_name).(region_name) additionally contains following
%  fields for all the pahses except baseline
%       - slope_mean: a matrix of size n_types x n_bands wich stores
%       the average slopes of that region and phase. the order of
%       the rows are sorted based on the visualization on
%       the figure, i.e. if exp is depicted first, the first row of the 
%       matrix is for exp data
%       - slope_se: a matrix of size n_types x n_bands wich stores
%       the standard error of the average slopes of that region and 
%       phase. the order of the rows are sorted based on the visualization 
%       on the figures, i.e. if exp is depicted first, the first row of the 
%       matrix is for exp data
%       - slope_p_val: a matrix of size n_comparasions x 1 which
%       sores the p-values of comparing slopes between differnt
%       data types. using un-paired ttest followed by FDR
% - spectrum_stats_out.compare.(type_name).(region_name) exists if the analyses was 
%   without pooling and contains the following fields for
%       region = region_str (i.e. PFC_rs,...) and type = type_name
%       (exp,ctrl,...):
%       - slope_mean: a matrix of size n_phases x n_bands wich stores
%       the average slopes of that region and data type for all phases. the order of
%       the rows are sorted based on the visualization on
%       the figure, i.e. if baseline is depicted first, the first row of the 
%       matrix is for exp data
%       - slope_se: a matrix of size n_phases x n_bands wich stores
%       the standard error of average slopes of that region and data type for all phases. the order of
%       the rows are sorted based on the visualization on
%       the figure, i.e. if baseline is depicted first, the first row of the 
%       matrix is for baseline data
%       - slope_p_val: a matrix of size n_comparasions x 1 which
%       sores the p-values of comparing slopes within same data type for differnt
%       phases. using paired ttest followed by FDR
%      - slope_individual: a cell array of length n_phases within each all
%      the individual slopes for all independent data of that data type and that region 
% -----------------------------------------------------------------------------
% this routine, loops over subjects, determines their type. concatenates
% the spectrum of independent samples from subjects of the same type. then 
% for each phase,region visualizes the avg +/- se of the %spectrum in log10 
% on a seperate figure for each data type. and overlays on top of them, the
% average estimated spectrum using the estimated line coefficients on the
% frequency range estimated. 
% then it takes the slope for each independt sample and compares them
% within each phase and each region between different dataypes using unpaired t-test followed by
% FDR correction. 
% Additionally, if there was no pooling over time was involved in the
% original analyses, it proceeds to compare the slopes for data of the same
% data type for each region across different phases using paired t-test
% followed by a FDR correction


%% initialize
n_subjects = length(results_all);
phases = fieldnames(results_all{1}.results.slope);
regions = fieldnames(results_all{1}.results.slope.(phases{1}));

n_phases = length(phases);
n_regions = length(regions);

type_all = cell(1,n_subjects);
for subj=1:n_subjects
    type_all{subj} = results_all{subj}.type;
end

types = unique(type_all,'stable');
n_types = length(types);
freq = slope.frequency;


%% plot spectrum
line_props = {'-b','-r','-g','-k'};
colors = {'blue','red','green','black'};

for phase=1:n_phases  
    figure('Name','Raw Spectrum with slope','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);sgtitle(phases{phase});
    for tp=1:n_types % for each data type
        for reg=1:n_regions % for each region            
            for t=1:n_types % initialize for each type
                tmp.(types{t}) = [];
                beta.(types{t}) = [];
            end
            
            for subj=1:n_subjects % loop over subjects
                type = results_all{subj}.type; % get data type
                tmp.(type) = horzcat(tmp.(type),results_all{subj}.results.slope.(phases{phase}).(regions{reg}).spectrum); % concatenate spectrum
                beta.(type) = horzcat(beta.(type),reshape(cell2mat(results_all{subj}.results.slope.(phases{phase}).(regions{reg}).slope),2,[])); % concatenate line parameters 
                f = results_all{subj}.results.slope.(phases{phase}).(regions{reg}).spectrum_freq; % get spectrum frequency
            end
            
            y = mean(log10(tmp.(types{tp})),2); % take the average spectrum in log10
            y_err = std(log10(tmp.(types{tp})),[],2)./sqrt(size(tmp.(types{tp}),2)); % get the S.E of spectrum average in log10
            x = f;
            idx_line_fit = f>freq(1) & f<freq(2); % find the index of the fitted line
            x_line = x(idx_line_fit); % get the frequencies of in the range of the fitted line
            y_line = polyval(mean(beta.(types{tp}),2),log10(x_line)); % evaluate the fitted line
            slope_stats_out.(phases{phase}).(regions{reg}).(types{tp}).spec_freq = f;
            slope_stats_out.(phases{phase}).(regions{reg}).(types{tp}).spec_mean = y;
            slope_stats_out.(phases{phase}).(regions{reg}).(types{tp}).spec_se = y_err;
            slope_stats_out.(phases{phase}).(regions{reg}).(types{tp}).beta = mean(beta.(types{tp}),2); 
            subplot(1,n_regions,reg);xlabel('frequency (Hz)');title(regions{reg})
            plot_shadedErrorBar(x,y,y_err,'lineProps',line_props{tp});
            hold all
            plot(x_line,y_line,'-m','LineWidth',3);
            axis square
            set_figure_properties()
        end
    end
    dim = [.9 .8 .2 .2];
    go_down = [0,0.035,0,0];
    for tp=1:n_types
        str = types{tp};
        dim_to_use = dim - (tp-1)*go_down;
        annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',colors{tp},'FontSize',12);
    end
    subplot(1,n_regions,1);ylabel('power');
    
end

%% comparing slopes
to_stats = struct;
for phase=1:n_phases % loop over phases  
    for reg=1:n_regions % loop over regions
        for t=1:n_types
            tmp.(types{t}) = []; % init data for each type
        end
        for subj=1:n_subjects % loop over subjects
            type = results_all{subj}.type; % get the type
            tmp.(type) = horzcat(tmp.(type),reshape(cell2mat(results_all{subj}.results.slope.(phases{phase}).(regions{reg}).slope),2,[])); % concatenate slopes for the data of the same type
        end
        to_stats.(phases{phase}).(regions{reg}) = tmp;
    end
    
end


n_comp = (n_types*(n_types -1))/2; 
for phase=1:n_phases % loop over phases
    figure('Name','Comparing Slopes','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);sgtitle(phases{phase})
    for reg=1:n_regions % loop over regions
        to_bar = zeros(n_types,1);
        err_bar = zeros(n_types,1);
        p_val = zeros(n_comp,1);
        
        to_test = cell(1,n_types);
        for type=1:n_types
            to_test{type} = to_stats.(phases{phase}).(regions{reg}).(types{type})(1,:); % get the data of that type
            to_bar(type,1) = mean(to_test{type}); % get the mean of slopes
            err_bar(type,1) = std(to_test{type})/sqrt(length(to_test{type})); % get the S.E of slopes
        end
        comp = 1;
        for type1=1:n_types % loop over type 1
            for type2=type1+1:n_types % loop over type 2
                [~,p_val(comp,1)] = ttest2(to_test{type1},to_test{type2}); % compate type 1 vs type 2
                comp = comp + 1;
            end
        end
        
        
        [~,~,~,p_val_corrected] = stat_fdr_bh(p_val); % correct p-value
        % plot results   
        subplot(1,n_regions,reg);ylim([-1,1]);
        bar(1:n_types,to_bar);hold all;
        errorbar(1:n_types,to_bar,err_bar,'.');
        xticklabels(types);
        comp = 1;
        for type1=1:n_types
            for type2=type1+1:n_types
                sigstar([type1,type2],p_val_corrected(comp));
            end
        end
        title(regions{reg});set_figure_properties()
        slope_stats_out.(phases{phase}).(regions{reg}).slope_mean = to_bar;
        slope_stats_out.(phases{phase}).(regions{reg}).slope_se = err_bar;
        slope_stats_out.(phases{phase}).(regions{reg}).slope_p_val = p_val_corrected;
    end
    subplot(1,n_regions,1);ylabel('fitted slope');
    
end

%% compare slopes within group
if results_all{1}.params.analyses.slope.pooling == 0 % only do within groups if there is no pooling
    n_comp = (n_phases*(n_phases -1))/2;
    for tp=1:n_types
        figure('Name','Comparing slope across conditions','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);sgtitle(types{tp});
        for reg=1:n_regions % for each region
            
            p_val = zeros(n_comp,1);
            to_bar = zeros(n_phases,1);
            err_bar = zeros(n_phases,1);
            for phase=1:n_phases
                to_bar(phase,:) = mean(to_stats.(phases{phase}).(regions{reg}).(types{tp})(1,:)); % get average slope for that phase
                err_bar(phase,:) = std(to_stats.(phases{phase}).(regions{reg}).(types{tp})(1,:))/sqrt(size(to_stats.(phases{phase}).(regions{reg}).(types{tp}),2)); % get S.E for that phase
            end
            comp = 1;
            for phase1=1:n_phases % for phase 1 
                for phase2=phase1+1:n_phases % for phase 2
                    to_test{phase1} = to_stats.(phases{phase1}).(regions{reg}).(types{tp})(1,:);
                    to_test{phase2} = to_stats.(phases{phase2}).(regions{reg}).(types{tp})(1,:);
                    [~,p_val(comp,1)] = ttest(to_test{phase1},to_test{phase2}); % compare slopes
                    
                    comp = comp + 1;
                end
            end
            
            [~,~,~,p_val_corrected] = stat_fdr_bh(p_val); % correct for p-val
            % plot results
            subplot(1,n_regions,reg);
            bar(1:n_phases,to_bar);hold all
            errorbar(1:n_phases,to_bar,err_bar,'.');
            comp = 1;
            for phase1=1:n_phases
                for phase2=phase1+1:n_phases
                    sigstar([phase1,phase2],p_val_corrected(comp));
                    comp = comp + 1;
                end
            end
            title(regions{reg});xticklabels(phases);set_figure_properties()
            slope_stats_out.compare.(types{tp}).(regions{reg}).slope_mean = to_bar;
            slope_stats_out.compare.(types{tp}).(regions{reg}).slope_se = err_bar;
            slope_stats_out.compare.(types{tp}).(regions{reg}).slope_p_val = p_val_corrected;
            slope_stats_out.compare.(types{tp}).(regions{reg}).slope_individual = to_test;
        end
        subplot(1,n_regions,1);ylabel('slope');
        
    end
end
end


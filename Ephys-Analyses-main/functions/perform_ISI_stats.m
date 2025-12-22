function [ISI_stats_out] = perform_ISI_stats(results_all)
% Function to do ISI stats
% Inputs: 
% Results_all: a cell array of length of n_subjects iside each cell there 
% is a structure with the following fields:
% - type: the type of the subject (exp, ctrl,...)
% - params: a structure containing the parameters of all the analyses. for
% example params.ISI contains the parameters of the ISI analyses
% - results: a structure containing the results of various analyses. for
% example results.ISI contains the results of the ISI analyses
% Outputs:
% ISI_stats_out: a structure containing various fields each with the
% results of the statistical analyses performed:
% - ISI_stats_out.(phase_name).(region_name), is a structure
% containing various info for the for phase = phase_name (baseline,cno,..) 
% and the region_name (PFC,...)  with the following fields:
%       - ISI_mean: a matrix of size n_types x 1 wich stores
%       the average ISI over all the channels of data of the same type of 
%       that region and phase. the order of
%       the rows are sorted based on the visualization on
%       the figure, i.e. if exp is depicted first, the first row of the 
%       matrix is for exp data
%       - ISI_se: a matrix of size n_types x 1 wich stores
%       the standard error of the average ISI over all the channels of data of the same type of 
%       that region and phase. the order of
%       the rows are sorted based on the visualization on
%       the figure, i.e. if exp is depicted first, the first row of the 
%       matrix is for exp data
%       - ISI_p_val: a matrix of size n_comparasions x 1 which
%       stores the p-values of comparing median of ISI between differnt
%       data types for that phase and region using ranksum test
%       - (type).histogram: ISI for all the individual channels of the data of the
%       type_name = type (exp,ctrl,...)
% - ISI_stats_out.phase_vs.(region_name).(type_name), is a structure
% containing various info for the region_name (PFC,...)  
% and type_name = type (exp,ctrl,...) with the following fields:
%       - ISI_mean: a matrix of size n_phases x 1 wich stores
%       the average ISI over all the channels of data of the same type of 
%       that region for each phase. the order of data is sorted based on order
%       of phases
%       - ISI_se: a matrix of size n_phases x 1 wich stores
%       the standard error of average ISI over all the channels of data of the same type of 
%       that region for each phase. the order of data is sorted based on order
%       of phases
%       - ISI_p_val: a matrix of size n_comparasions x 1 which
%       stores the p-values of comparing median of ISI between differnt
%      phases for that datatype and region using ranksum test
% -----------------------------------------------------------------------------
% this routine, will loop ovar subjects and concatenates the calculated
% ISI values for all channels of subects of the same type together. then it starts by 
% visualizing the histogram of ISIs for each data type, region and phase.
% Then it proceeds to first compare ISI of different groups for each phase
% and region using ranksum followed by FDR and then compare the ISI between phases within
% same groups for each region using ranksum followed by FDR. 
%% initialize
n_subjects = length(results_all);
phases = fieldnames(results_all{1}.MUA_results.results.rate);
regions = fieldnames(results_all{1}.MUA_results.results.rate.(phases{1}));

n_phases = length(phases);
n_regions = length(regions);

type_all = cell(1,n_subjects);
for subj=1:n_subjects
    type_all{subj} = results_all{subj}.MUA_results.type;
end

types = unique(type_all,'stable');
n_types = length(types);

ISI_stats_out = struct;

%% plot ISI histogram
to_stats = struct;
edges = 0:0.002:2;
for tp=1:n_types % loop over data types
    figure('Name','Inter Spike Interval Histogram','NumberTitle','off','units','normalized','outerposition',[0 0 95 0.95]);sgtitle(types{tp});
    count = 1;
    for phase=1:n_phases % loop over phases
        for reg=1:n_regions 
            for t=1:n_types
                tmp.(types{t}) = [];
            end
            for subj=1:n_subjects
                type = results_all{subj}.MUA_results.type; % get data type
                tmp.(type) = horzcat(tmp.(type),results_all{subj}.MUA_results.results.ISI.(phases{phase}).(regions{reg}).ISI); % concatenate data of same type
            end
        ISI = tmp.(types{tp});
        n_units = length(ISI);
        to_hist = [];
        for unit=1:n_units % loop over channels
            to_stats.(phases{phase}).(regions{reg}).(types{tp})(unit) = mean(ISI{unit}); % take the average of channel
            to_hist = [to_hist;ISI{unit}];
        end
        subplot(n_phases,n_regions,count);
        count = count + 1;
        histogram(to_hist,edges,'Normalization','probability');xline(mean(to_hist),'--r','linewidth',2);title(strcat(phases{phase},'-',regions{reg}));
        ISI_stats_out.(phases{phase}).(regions{reg}).(types{tp}).histogram = to_hist;
        ylim([0,0.05]);
        axis square
        set_figure_properties()
        end
    end
end


%% within group across phase statistics
n_comp = (n_phases*(n_phases -1))/2;
for type=1:n_types % loop over data type
     figure('Name','Comparing ISI within Groups','NumberTitle','off','units','normalized','outerposition',[0 0 0.95 0.95]);sgtitle(types{type});
    for reg=1:n_regions % loop over regions 
        to_bar = zeros(n_phases,1);
        err_bar = zeros(n_phases,1);
        p_val = zeros(n_comp,1);
        to_test = cell(1,n_phases);
        for phase=1:n_phases
            to_test{phase} = to_stats.(phases{phase}).(regions{reg}).(types{type});
            to_bar(phase,1) = mean(to_test{phase},'omitnan'); % get the average
            err_bar(phase,1) = std(to_test{phase},'omitnan')/sqrt(length(to_test{type})); % get the S.E
        end
        
        
        comp = 1;
        for phase1=1:n_phases % loop over phase1
            for phase2=phase1+1:n_phases % loop over phase2
                p_val(comp,1) = ranksum(to_test{phase1},to_test{phase2}); % compare
                comp = comp + 1;
            end
        end
        
        [~,~,~,p_val_corrected] = stat_fdr_bh(p_val); % correct p-val
        % plot
        subplot(1,n_regions,reg);
        bar(1:n_phases,to_bar);hold all;
        errorbar(1:n_phases,to_bar,err_bar,'.');
        ylim([0,2]);xticklabels(phases);
        set_figure_properties()
        comp = 1;
        for phase1=1:n_phases
            for phase2=phase1+1:n_types
                sigstar([phase1,phase2],p_val_corrected(comp));
                comp = comp+1;
            end
        end
        title(regions{reg});
        ISI_stats_out.phase_vs.(regions{reg}).(types{type}).ISI_mean = to_bar;
        ISI_stats_out.phase_vs.(regions{reg}).(types{type}).ISI_se = err_bar;
        ISI_stats_out.phase_vs.(regions{reg}).(types{type}).ISI_p_val = p_val_corrected;
        
        subplot(1,n_regions,1);ylabel('Inter Spike Interval mean');       
    end
end


%% across group within same phase statistics
n_comp = (n_types*(n_types -1))/2;
for phase=1:n_phases % loop over phases
    figure('Name','Comparing ISI across groups','NumberTitle','off','units','normalized','outerposition',[0 0 95 0.95]);sgtitle(phases{phase})
    for reg=1:n_regions % loop over region
        to_bar = zeros(n_types,1);
        err_bar = zeros(n_types,1);
        p_val = zeros(n_comp,1);
        to_test = cell(1,n_types);
        for type=1:n_types % loop over data types
            to_test{type} = to_stats.(phases{phase}).(regions{reg}).(types{type}); % get the data of that group
            to_bar(type,1) = mean(to_test{type},'omitnan'); % get the average
            err_bar(type,1) = std(to_test{type},'omitnan')/sqrt(length(to_test{type})); % get the SE
        end
        comp = 1;
        for type1=1:n_types % loop type 1
            for type2=type1+1:n_types % loop type 2
                p_val(comp,1) = ranksum(to_test{type1},to_test{type2}); % compare across types
                comp = comp + 1;
            end
        end
        
        
        [~,~,~,p_val_corrected] = stat_fdr_bh(p_val); % correct p-val
        % plot    
        subplot(1,n_regions,reg);ylim([-1,1]);
        bar(1:n_types,to_bar);hold all;
        errorbar(1:n_types,to_bar,err_bar,'.');
        ylim([-1,1])
        xticklabels(types);
        set_figure_properties()
        comp = 1;
        for type1=1:n_types
            for type2=type1+1:n_types
                sigstar([type1,type2],p_val_corrected(comp));
                comp = comp + 1;
            end
        end
        title(regions{reg});
        ISI_stats_out.(phases{phase}).(regions{reg}).ISI_mean = to_bar;
        ISI_stats_out.(phases{phase}).(regions{reg}).ISI_se = err_bar;
        ISI_stats_out.(phases{phase}).(regions{reg}).ISI_p_val = p_val_corrected;
    end
    subplot(1,n_regions,1);ylabel('mean Inter Spike Interval');
    
end
end


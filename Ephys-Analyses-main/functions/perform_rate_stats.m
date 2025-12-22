function rate_stats_out = perform_rate_stats(results_all,rt)
% Function to do rate stats
% Inputs: 
% Results_all: a cell array of length of n_subjects iside each cell there 
% is a structure with the following fields:
% - type: the type of the subject (exp, ctrl,...)
% - params: a structure containing the parameters of all the analyses. for
% example params.rate contains the parameters of the rate analyses
% - results: a structure containing the results of various analyses. for
% example results.rate contains the results of the rate analyses
% rt: a structure containing some details for rate analyses with the
% following fields:
%  - ch_lvl: 1 to do ch_lvl analyses, 0 for channel average 
%  - pooling an integer value to pool data in minutes, 0 for no pooling
%  - eff: 'lower' if u expect rate going low, 'higher' if u expect it to go high
% Outputs:
% rate_stats_out: a structure containing various fields each with the
% results of the statistical analyses performed:
% - rate_stats_out.(phase_name).(region_name), is a structure
% containing various info for the for phase = phase_name except baseline
% (cno,..) and the region_name (PFC,...)  with the following fields:
%       - modulation_mean: a matrix of size n_types x 1 wich stores
%       the average rate modulation index of that region and phase. the order of
%       the rows are sorted based on the visualization on
%       the figure, i.e. if exp is depicted first, the first row of the 
%       matrix is for exp data
%       - modulation_se: a matrix of size n_types x 1 wich stores
%       the standard error of the average rate modulation index of that region and 
%       phase. the order of the rows are sorted based on the visualization 
%       on the figures, i.e. if exp is depicted first, the first row of the 
%       matrix is for exp data
%       - modulation_p_val: a matrix of size n_comparasions x 1 which
%       sores the p-values of comparing modulation indices between differnt
%       data types for that phase and region
%       - modulation_individual:  a cell array n_types x 1  (data types
%       (exp,ctrl,...))
%       inside each cell all the independet rate modulations
%       of that phase,region and data type is stored in columns of a
%       matrix. Note that since pooiling was involved pooled data of a
%       subjects are considered independent. Also, the first n columns 
%       are the first n pooled segments of a subject, then the next subject
%       the order of the data in the cell are sorted based on the visualization 
%       on the figures, i.e. if exp is depicted first, the first element of cell is for exp data
% - rate_stats_out.(region_name).(type_name), is a structure
% containing various info for the region_name (PFC,...)  
% and type_name = type (exp,ctrl,...) with the following fields:
%       - ch_avg_rate_mean: an array of length n_timepoints which the
%       grand average of rate (across channels and subjects) for each
%       timepoint is stored. 
%       - ch_avg_rate_se: an array of length n_timepoints inside the
%       standard error of grand average of rate (across channels and subjects) for each
%       timepoint is stored. 
% -----------------------------------------------------------------------------
% this routine, will loop ovar subjects and concatenates the calculated
% spike rate values for data of the same type together. then it starts by 
% visualizing the grand average (average over subject and channel) of spike rate 
% for each datatype (exp,ctrl) over time on a seperate figure. Additionally
% if there are only 2 phases, it will scatter plot the average over time of rate in phase2
% vs phase1 for all the channels of data of the same type seperately for each region.
% then it checks if pooling over channel is requested, if not it takes the
% average over channels for each subject at each time point. otherwise it
% will consider the channels as independent. if pooling over time is
% requested, it will also pool the data resulted from the previous part
% over requested time points and treats them as independent. 
% it proceeds by calculating a modulation index for all the independent samples for all
% phases as follow: 
% MI = (rate(phase_2,region_x)/average(rate(baseline,region)),average over
% time))
% then checks whether the median of modulation index differens across data
% types for each phase (except baseline) and region using ranksum test
% followed by a FDR correction.
% Note: it will always treat the phase1 as baseline
% Note: it will plot the grand average of rate as if the data is continous
% over time whereas in reality the timing might not be continous, due to
% requests in the original analyses. 
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
if ~strcmp(types{1},'ctrl')
    types{2} = types{1};
    types{1} = 'ctrl';
end

rate_stats_out = struct;
%% plot spike rate
line_props = {'-b','-r','-g','-k'}; % colors for plotting 
tycolors.ctrl = 'blue';
tycolors.exp = 'red';
phcolors.baseline = '#1AFF1A';
phcolors.cno=  '#FF00FF';
phcolors.CNO=  '#FF00FF';

figure('Name','Raw spike rate','NumberTitle','off','units','normalized','outerposition',[0 0 9 0.9]);
st = sgtitle("Raw spike rate");
st.FontSize = 30;
for tp=1:n_types % loop over type
    %figure('Name','Raw spike rate','NumberTitle','off','units','normalized','outerposition',[0 0 9 0.9]);
    %st = sgtitle(types{tp});
    %st.FontSize = 30;
    for reg=1:n_regions % loop over region
        for phase=1:n_phases       
            for t=1:n_types
                tmp.(types{t}).all = [];
                tmp.(types{t}).indmean = [];
            end
            for subj=1:n_subjects % loop over subjects
                type = results_all{subj}.MUA_results.type; % get data type
                tmp.(type).all = horzcat(tmp.(type).all,results_all{subj}.MUA_results.results.rate.(phases{phase}).(regions{reg}).rate); % concatenate data of the same type
                tmp.(type).indmean = horzcat(tmp.(type).indmean,mean(results_all{subj}.MUA_results.results.rate.(phases{phase}).(regions{reg}).rate,2));
            end
            if phase == 1
               x = []; 
               y = [];
               y_err = [];
               indmean = [];
            else 
                t_fill = [1:results_all{1}.MUA_results.time.(phases{phase}).ti] +length(x);
                x = [x,t_fill];               
                y = [y;zeros(results_all{1}.MUA_results.time.(phases{phase}).ti,1)];
                y_err = [y_err;zeros(results_all{1}.MUA_results.time.(phases{phase}).ti,1)];
            end
            y = [y;mean(tmp.(types{tp}).all,2)]; % get the average and concate
            y_err = [y_err;std(tmp.(types{tp}).all,[],2)./sqrt(size(tmp.(types{tp}).all,2))];
            x = [x,[1:size(tmp.(types{tp}).all,1)] + (phase - 1)*length(x)];
            indmean = vertcat(indmean,tmp.(types{tp}).indmean);
            subplot(1,n_regions,reg);title(regions{reg});
            xline(x(end),'Color',phcolors.(phases{phase}),'linewidth',2,'Label',phases{phase});
            xlabel('time');
        end
        subplot(1,n_regions,reg);title(regions{reg});
        plot_shadedErrorBar(x,y,y_err,'lineprops',tycolors.(types{tp}));xlim([x(1),x(end)+1]);
        rate_stats_out.(regions{reg}).(types{tp}).ch_avg_rate_mean = y;
        rate_stats_out.(regions{reg}).(types{tp}).ch_avg_rate_se = y_err;
        rate_stats_out.(regions{reg}).(types{tp}).indmean = indmean;
        axis square
        set_figure_properties()
    end
    dim = [.9 .8 .2 .2];
    go_down = [0,0.035,0,0];
    for type=1:n_types
        str = types{type};
        dim_to_use = dim - (type-1)*go_down;
        annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',tycolors.(types{type}),'FontSize',20);
    end
end


%% plotting the phase1 vs phase2 for all channels 
% average over time for each subject and each channel
if n_phases == 2 % only do if there are 2 phases
    for tp=1:n_types
        figure('Name','pre-vs-post','NumberTitle','off','units','normalized','outerposition',[0 0 9 0.9]);
        st = sgtitle(types{tp});
        st.FontSize = 30;
        for reg=1:n_regions
            to_scatter = [];
            for phase=1:n_phases
                for t=1:n_types
                    tmp.(types{t}).all = [];
                    tmp.(types{t}).indmean = [];
                end
                for subj=1:n_subjects
                    type = results_all{subj}.MUA_results.type;
                    tmp.(type).all = horzcat(tmp.(type).all,results_all{subj}.MUA_results.results.rate.(phases{phase}).(regions{reg}).rate);
                end
                to_scatter(:,phase) = mean(tmp.(types{tp}).all);     
            end
            subplot(1,n_regions,reg);
            scatter(to_scatter(:,1),to_scatter(:,2))
            max_rate = max(to_scatter(:));
            axis([0 max_rate+1 0 max_rate+1])
            hold all 
            plot(0:0.1:max_rate+1,0:0.1:max_rate+1,'--r')
            if strcmp(rt.eff,'lower')
                per_low = (sum(to_scatter(:,2) < to_scatter(:,1))./size(to_scatter,1))*100;
                str_ant = strcat(num2str(per_low),'%','decreasing');
            elseif strcmp(rt.eff,'higher')
                per_high = (sum(to_scatter(:,2) > to_scatter(:,1))./size(to_scatter,1))*100;
                str_ant = strcat(num2str(per_high),' %',' increasing');
            end           
            xlabel(phases{1});ylabel(phases{2})
            title(regions{reg});
            set_figure_properties()
            x_pos = 0.2 + (reg-1)*0.5;
            annotation('textbox',[x_pos 0.5 0.3 0.3],'String',str_ant,'FitBoxToText','on');
        end
    end
    
end

%Compute individual rates
to_stats = struct;
for phase=2:n_phases
    figure('Name','Spike rate single subject','NumberTitle','off','units','normalized','outerposition',[0 0 0.95 0.95]);
    st = sgtitle(phases{phase});
    st.FontSize = 30;
    for reg=1:n_regions
        for t=1:n_types
            tmp.(types{t}).all = [];
            tmp.(types{t}).indmean = [];
        end
        for subj=1:n_subjects
            type = results_all{subj}.MUA_results.type;
            if rt.ch_lvl == 0
                tmp_rate = mean(results_all{subj}.MUA_results.results.rate.(phases{phase}).(regions{reg}).rate,2);
                tmp_baseline = mean(mean(results_all{subj}.MUA_results.results.rate.(phases{1}).(regions{reg}).rate,2),1);
            else
                tmp_rate = results_all{subj}.MUA_results.results.rate.(phases{phase}).(regions{reg}).rate;
                tmp_baseline = mean(results_all{subj}.MUA_results.results.rate.(phases{1}).(regions{reg}).rate,1);
            end
            
            if rt.pooling >0
                tmp_rate_pool = reshape(tmp_rate,[],rt.pooling,size(tmp_rate,2));
                tmp_rate_pool = squeeze(mean(tmp_rate_pool,2));
            else
                tmp_rate_pool = mean(tmp_rate);
            end
            tmp_subj_normalized = tmp_rate_pool./(tmp_baseline);
            tmp.(type).all = cat(2,tmp.(type).all,tmp_subj_normalized);
         
            %Plot single subject firing rates (normalized to baseline)
            subplot(1,n_regions,reg);title(regions{reg});
            ti = results_all{1}.MUA_results.time.(phases{phase}).ti;
            x = ti:(ti+size(tmp_rate,1)-1);
            y = tmp_rate./(tmp_baseline);
            %tmp.(type).indmean = cat(2,tmp.(type).indmean,y);
            y_err = std(y,[],2)/sqrt(size(y,2));
            ym = mean(y,2);
            tmp.(type).indmean = cat(2,tmp.(type).indmean,y);
            plot(x,ym,'Color',tycolors.(type)); hold all;
        end
        xlabel('time (mins after injection)');
        ylabel('baseline normalized rate')
        set_figure_properties();
        to_stats.(phases{phase}).(regions{reg}) = tmp;
        for t=1:n_types
            rate_stats_out.(regions{reg}).(types{t}).(phases{phase}).rate_normalized = tmp.(types{t}).indmean;
        end
    end
    dim = [.9 .8 .2 .2];
    go_down = [0,0.035,0,0];
    for type=1:n_types
        str = types{type};
        dim_to_use = dim - (type-1)*go_down;
        annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',tycolors.(types{type}),'FontSize',20);
    end
end
%Plot baseline normalized firing rate: 
figure('Name','Baseline normalized firng rate','NumberTitle','off','units','normalized','outerposition',[0 0 0.95 0.95]);
y = mean(rate_stats_out.PFC.ctrl.cno.rate_normalized,2);
y_err = std(rate_stats_out.PFC.ctrl.cno.rate_normalized,[],2)./sqrt(size(rate_stats_out.PFC.ctrl.cno.rate_normalized,2));
plot_shadedErrorBar(x,y,y_err,'lineprops','-b');
y = mean(rate_stats_out.PFC.exp.cno.rate_normalized,2);
y_err = std(rate_stats_out.PFC.exp.cno.rate_normalized,[],2)/sqrt(size(rate_stats_out.PFC.exp.cno.rate_normalized,2));
hold on
plot_shadedErrorBar(x,y,y_err,'lineprops','-r');
xlabel('time (mins after injection)');
ylabel('baseline normalized rate')
dim = [.9 .8 .2 .2];
go_down = [0,0.035,0,0];
for type=1:n_types
str = types{type};
dim_to_use = dim - (type-1)*go_down;
annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',tycolors.(types{type}),'FontSize',20);
end
set_figure_properties()

%% compute the modulation index
n_comp = (n_types*(n_types -1))/2;
for phase=2:n_phases
    figure('Name','Spike rate Modulation','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    st = sgtitle(phases{phase});
    st.FontSize = 30;
    for reg=1:n_regions
        to_bar = zeros(n_types,1);
        err_bar = zeros(n_types,1);
        p_val = zeros(n_comp,1);
        
        to_test = cell(1,n_types);
        for type=1:n_types
            to_test{type} = to_stats.(phases{phase}).(regions{reg}).(types{type}).all(:);
            to_bar(type,1) = mean(to_test{type}(:));
            err_bar(type,1) = std(to_test{type}(:))/sqrt(length(to_test{type}(:)));
            rate_stats_out.(regions{reg}).(types{type}).(phases{phase}).modulation_individual = to_test{type};
        end
        comp = 1;
        for type1=1:n_types
            for type2=type1+1:n_types
                p_val(comp,1) = ranksum(to_test{type1},to_test{type2});
                comp = comp + 1;
            end
        end
        
        
        [~,~,~,p_val_corrected] = stat_fdr_bh(p_val);
           
        subplot(1,n_regions,reg);ylim([0,2]);
        bar(1:n_types,to_bar);hold all;
        errorbar(1:n_types,to_bar,err_bar,'.');
        ylim([0,2])
        xticklabels(types);       
        comp = 1;
        for type1=1:n_types
            for type2=type1+1:n_types
                sigstar([type1,type2],p_val_corrected(comp));
                comp = comp + 1;
            end
            rate_stats_out.(regions{reg}).(types{type1}).(phases{phase}).modulation_mean = to_bar(type1);
            rate_stats_out.(regions{reg}).(types{type1}).(phases{phase}).modulation_se = err_bar(type1);
        end
        title(regions{reg});set_figure_properties()
        rate_stats_out.(regions{reg}).modulation_across_types.(phases{phase}).modulation_p_val = p_val_corrected;
    end
    subplot(1,n_regions,1);ylabel('baseline normalized rate index');
end

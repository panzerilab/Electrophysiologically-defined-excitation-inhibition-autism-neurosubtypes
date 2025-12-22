function PLV_stats_out = perform_PLV_MUA_stats(results_all)
% Function to do PLV stats
% Inputs: 
% Results_all: a cell array of length of n_subjects iside each cell there 
% is a structure with the following fields:
% - type: the type of the subject (exp, ctrl,...)
% - params: a structure containing the parameters of all the analyses. for
% example params.ISI contains the parameters of the ISI analyses
% - results: a structure containing the results of various analyses. for
% example results.PLV contains the results of the PLV analyses
% Outputs:
% PLV_stats_out: a structure containing various fields each with the
% results of the statistical analyses performed:
%  - PLV_stats_out.compare.(phase_name).(region_name) contains following
%       fields for the phase = phase_name (except the baseline) and the
%       region = region_name (PFC,rs,...)
%       - modulation_mean: a matrix of size n_types x 1 wich stores
%       the average modulation index of that region and phase. the order of
%       the rows are sorted based on the visualization on
%       the figure, i.e. if exp is depicted first, the first row of the 
%       matrix is for exp data
%       - modulation_se: a matrix of size n_types x 1 wich stores
%       the standard error of the average modulation index of that region and 
%       phase. the order of the rows are sorted based on the visualization 
%       on the figures, i.e. if exp is depicted first, the first row of the 
%       matrix is for exp data
%       - modulation_p_val: a matrix of size n_comparasions x 1 which
%       sores the p-values of comparing modulation indices between differnt
%       data types.
%   - PLV_stats_out.(type_name).(phase_name).(region_name).plv_individual: 
%     contains the individual PLV for all the channels with plv>0.1 for the
%     data of type = type_name (exp,ctrl,...), phase = phase_name
%     (baseline,cno,...) and region = region_name (PFC,rs,...)
% -----------------------------------------------------------------------------
% This function loops over different datatypes, first visualizes the preffered angles 
% for the data with the plv>0.1 on a cosine seperately for each data type and region and phase.
% then plots the histogram of this preferred phases. Then proceeds to
% calculate a modulation index as follow:
% note: it will always consider the first phase as baseline
% MI = (PLV(region_x,phase_y) -
% PLV(region_x,baseline)/(PLV(region_x,phase_y) + PLV(region_x,baseline)
% then moves to compare the PLV of different phases and regions across data
% types using ranksum followed by FDR
% Requires circ_stat toolbox   
%% initialize
n_subjects = length(results_all);
phases = fieldnames(results_all{1}.MUA_results.results.PLV);
regions = fieldnames(results_all{1}.MUA_results.results.PLV.(phases{1}));

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
n_subjects = length(results_all);
n_phases = length(phases);
n_regions = length(regions);




%% visualizing the preferred phase
x = 0:0.01:2*pi; 
for tp=1:n_types % loop over data types
    figure('Name','Units Preferred Phase of firing','NumberTitle','off','units','normalized','outerposition',[0 0 95 95]);
    st = sgtitle(types{tp});
    st.FontSize = 30;
    count = 1;
    for phase=1:n_phases % loop over phase
        for reg=1:n_regions  % loop over region 
            for t=1:n_types % init data of type
                tmp.(types{t}) = [];
                tmp2.(types{t}) = [];
            end
            for subj=1:n_subjects % loop over subject
                type = results_all{subj}.MUA_results.type; % get subject type
                tmp.(type) = horzcat(tmp.(type),results_all{subj}.MUA_results.results.PLV.(phases{phase}).(regions{reg}).mean); % cocatenate data of same type
                tmp2.(type) = horzcat(tmp2.(type),results_all{subj}.MUA_results.results.PLV.(phases{phase}).(regions{reg}).PLV); % cocatenate data of same type
                str = sprintf('%.0f,',results_all{subj}.MUA_results.results.PLV.(phases{phase}).(regions{reg}).ch_to_flip);
                disp(['Channels to flip in phase ',phases{phase},' in region ',regions{reg},' for subj ', num2str(subj) ,':',str ])
                
            end
        to_plot_mean = wrapTo2Pi(tmp.(types{tp})); % wrap data to 2pi
        to_plot_PLV = tmp2.(types{tp});
        to_plot = to_plot_mean(to_plot_PLV>0.1); % check threshold
        subplot(n_phases,n_regions,count);
        count = count + 1;
        edges = 0:pi/16:2*pi;polarhistogram(to_plot,edges,'normalization','probability');
        %plot(x,cos(x),'LineWidth',2);hold all;plot(to_plot,cos(to_plot)+0.1,'x');
        %axis square;
        title(strcat(phases{phase},'-',regions{reg}));
        set_figure_properties()
        PLV_stats_out.(types{tp}).(phases{phase}).(regions{reg}).PLV_individual = to_plot_PLV;
        end
    end
end

%% histogram of phases
edges = -pi/8:pi/4:2*pi+pi/8; % limits and edges of histogram
centers = movmean(edges,2); % center of bins of histogram
centers = centers(2:end-1); 
for tp=1:n_types % loop over type
    figure('Name','Histogram of phases of firing','NumberTitle','off','units','normalized','outerposition',[0 0 95 95]);
    st = sgtitle(types{tp});
    st.FontSize = 30;
    count = 1;
    for phase=1:n_phases % loop over phase
        for reg=1:n_regions % loop over region
            for t=1:n_types
                tmp.(types{t}) = [];
            end
            for subj=1:n_subjects % loop over subject
                type = results_all{subj}.MUA_results.type;
                tmp.(type) = horzcat(tmp.(type),results_all{subj}.MUA_results.results.PLV.(phases{phase}).(regions{reg}).phases); % cocatenate data of same type
                
            end
        to_hist = wrapTo2Pi(cell2mat(tmp.(types{tp})(:))); % wrap to 2pi
        [counts,edges] = histcounts(to_hist,edges); % get the counts
        counts(1) = counts(1)+counts(end);
        counts(end) = [];
        prob = counts./sum(counts); % calculate probability
        subplot(n_phases,n_regions,count);
        count = count + 1;
        bar(rad2deg(centers),prob);xticks(rad2deg(centers));
        [alpha,kappa] = circ_vmpar(to_hist); 
        PLV_stats_out.(types{tp}).(phases{phase}).(regions{reg}).vmpars = [alpha,kappa];
        %edges = 0:pi/16:2*pi;polarhistogram(to_hist,edges,'normalization','probability');
        axis square;
        title(strcat(phases{phase},'-',regions{reg}),['\alpha= ',num2str(alpha, '%.2f'),' \kappa = ',num2str(kappa, '%.2f')]);
        set_figure_properties()
        PLV_stats_out.(types{tp}).(phases{phase}).(regions{reg}).all_phases = to_hist;

        end
        
    end
end

%%Circular variance test
%Compare the distribution of phases across the two phases
% for tp = 1:n_types
%     for reg =1:n_regions
%         [pval, f] = circ_ktest(PLV_stats_out.(types{tp}).(phases{1}).(regions{reg}).all_phases, PLV_stats_out.(types{tp}).(phases{2}).(regions{reg}).all_phases);
%         PLV_stats_out.(types{tp}).(regions{reg}).kpval = pval;
%     end
% end

%% calculating normalized PLV and statistics

to_stats = struct;
for phase=2:n_phases % loop over phases except baseline
    for reg=1:n_regions % loop over regions
        for t=1:n_types
            tmp.(types{t}) = [];
        end
        for subj=1:n_subjects % loop over subjects
            type = results_all{subj}.MUA_results.type; % get data type
            tmp_PLV = results_all{subj}.MUA_results.results.PLV.(phases{phase}).(regions{reg}).PLV'; % get the data of that phase
            tmp_baseline = results_all{subj}.MUA_results.results.PLV.(phases{1}).(regions{reg}).PLV'; % get the baseline
            tmp_subj_normalized = (tmp_PLV - tmp_baseline)./((tmp_PLV + tmp_baseline + eps)); % calculate modulation index
            tmp.(type) = cat(1,tmp.(type),tmp_subj_normalized);
        end
        to_stats.(phases{phase}).(regions{reg}) = tmp;
    end
    
end


n_comp = (n_types*(n_types -1))/2;
for phase=2:n_phases % loop over phases except baseline
    figure('Name','Spike-LFP PLV Modulation','NumberTitle','off','units','normalized','outerposition',[0 0 0.95 0.95]);
    st = sgtitle(phases{phase});
    st.FontSize = 30;
    for reg=1:n_regions % loop over regions
        to_bar = zeros(n_types,1);
        err_bar = zeros(n_types,1);
        p_val = zeros(n_comp,1);
        to_test = cell(1,n_types);
        for type=1:n_types % loop over types
            to_test{type} = to_stats.(phases{phase}).(regions{reg}).(types{type}); 
            to_bar(type,1) = mean(to_test{type},'omitnan'); % get the mean
            err_bar(type,1) = std(to_test{type},'omitnan')/sqrt(length(to_test{type})); % get the S.E
        end
        comp = 1;
        for type1=1:n_types % loop over type 1 
            for type2=type1+1:n_types % loop over type 2
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
        comp = 1;
        for type1=1:n_types
            for type2=type1+1:n_types
                sigstar([type1,type2],p_val_corrected(comp));
                comp = comp + 1;
            end
        end
        title(regions{reg});set_figure_properties()
        PLV_stats_out.compare.(phases{phase}).(regions{reg}).modulation_individual = to_test; 
        PLV_stats_out.compare.(phases{phase}).(regions{reg}).modulation_mean = to_bar;
        PLV_stats_out.compare.(phases{phase}).(regions{reg}).modulation_se = err_bar;
        PLV_stats_out.compare.(phases{phase}).(regions{reg}).modulation_p_val = p_val_corrected;
    end
    subplot(1,n_regions,1);ylabel('PLV modulation index');
end


end


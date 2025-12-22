function BLP_stats_out = perform_BLP_stats(results_all)
% Function to do BLP stats
% Inputs: 
% Results_all: a cell array of length of n_subjects iside each cell there 
% is a structure with the following fields:
% - type: the type of the subject (exp, ctrl,...)
% - params: a structure containing the parameters of all the analyses. for
% example params.BLP contains the parameters of the BLP analyses
% - results: a structure containing the results of various analyses. for
% example results.BLP contains the results of the BLP analyses
% Outputs:
% BLP_stats_out: a structure containing various fields each with the
% results of the statistical analyses performed:
% - BLP_stats_out.(phase_name).(region_str), is a structure
% containing various info for the for phase = phase_name
% (baseline,cno,..) and the region_str (PFC-rs,...)  with the following fields:
%   - BLP_mean: an array of size n_types x 1 wich stores
%       the the average BLP correlation of that region and 
%       phase. the order of the rows are sorted based on the visualization 
%       on the figures, i.e. if exp is depicted first, the first row of the 
%       matrix is for exp data
%    - BLP_se: an array of size n_types x 1 wich stores
%       the standard error of the average BLP correlation of that region and 
%       phase. the order of the rows are sorted based on the visualization 
%       on the figures, i.e. if exp is depicted first, the first row of the 
%       matrix is for exp data
%    - BLP_p_val: an array of size n_comparasions x 1 which
%       sores the p-values of comparing BLP correlation value between differnt
%       data types.
% - BLP_stats_out.corr_vs.(region_str).(type_name) exists if the analyses was 
%   without pooling and contains the following fields for
%       region = region_str (i.e. PFC_rs,...) and type = type_name
%       (exp,ctrl,...):
%       - mean: an array the length of n_phases which contains the
%        average of BLP correlation values for specific data type in the
%       region_str for each phase
%       - se: an array the length of n_phases which contains the
%        standard error of BLP correlation values for specific data type in the
%       region_str for each phase
%       - pval: an array containing the p-values for statistical
%       comparisions of change in average of BLP averages within the same
%       data type and region between different phases
% - BLP_stats_out.compare.(phase_name).(region_str) exists if the analyses was 
%   with pooling contains following fields for the phase = phase_name (except the baseline) and the
%       region = region_str (PFC_rs) 
%       - modulation_mean: a matrix of size n_types x 1 wich stores
%       the average BLP modulation index of that region and phase. the order of
%       the rows are sorted based on the visualization on
%       the figure, i.e. if exp is depicted first, the first row of the 
%       matrix is for exp data
%       - modulation_se: a matrix of size n_types x 1 wich stores
%       the standard error of the average BLP modulation index of that region and 
%       phase. the order of the rows are sorted based on the visualization 
%       on the figures, i.e. if exp is depicted first, the first row of the 
%       matrix is for exp data
%       - modulation_p_val: a matrix of size n_comparasions x 1 which
%       sores the p-values of comparing modulation indices between differnt
%       data types.
%       - modulation_individual:  a structure wih fieldnames corresponding
%       to data types (exp,ctrl,...) inside each all the independet BLP modulations
%       of that phase,region and data type is stored in columns of a
%       matrix. Note that since pooiling was involved pooled data of a
%       subjects are considered independent. Also, the first n columns 
%       are the first n pooled segments of a subject, then the next subject
%       
% -----------------------------------------------------------------------------
% this routine, will loop ovar subjects and concatenates the BLP
% correlation values for data of the same type together. Then if pooling
% was not requested in the analyses, it checks to see whether the average of
% BLP is different across phases within the same data type (exp,ctrl,...)
% and region (PFC-rs,...) using paired t-test followed by FDR correction. 
% if pooling was requested in the original analyses, it proceeds by
% calculating a modulation index for all the independent samples for all
% phases as follow: 
% MI = (BLP(phase_2,region_xy) - BLP(baseline,region_xy))/ (BLP(phase_2,region_xy) + BLP(baseline,region_xy))
% then checks whether the median of modulation index differens across data
% types for each phase (except baseline) and region using ranksum test
% followed by a FDR correction. Additionally, it tests whether the actual
% average (not moudlation index) varies across different data types for
% each phase and region using un-paired ttest followed by FDR.
% note: it will always consider the first phase as baseline
% note: BLP is always measure between two regions, so by region we
% mean region-pair (for example PFC-rs)

%% initialize
n_subjects = length(results_all);
phases = fieldnames(results_all{1}.results.BLP);
regions = fieldnames(results_all{1}.results.BLP.(phases{1}));
type_all = cell(1,n_subjects);
for subj=1:n_subjects
    type_all{subj} = results_all{subj}.type;
end
types = unique(type_all,'stable'); % get data type
n_types = length(types);
%Make sure first type is always ctrl
if ~strcmp(types{1},'ctrl')
    types{2} = types{1};
    types{1} = 'ctrl';
end

%% plotting coherency
line_props = {'-b','-r','-g','-k'}; % colors for plotting 
tycolors.ctrl = "#0072BD";
tycolors.exp = "#D95319";
phcolors.baseline = '#1AFF1A';
phcolors.cno=  '#FF00FF';
phcolors.CNO=  '#FF00FF';


n_subjects = length(results_all);
n_phases = length(phases);
n_regions = length(regions);
% fix region names from reg1_reg2 to reg1-reg2 for titles
regions_name = regions;
for reg=1:n_regions
    temp_reg = regions{reg};
    idx_underline = strfind(temp_reg,'_');
    temp_reg(idx_underline) = '-';
    regions_name{reg} = temp_reg;
end




for tp=1:n_types % loop over types
    figure('Name','Raw BLP','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    st = sgtitle(types{tp});
    st.FontSize = 30;
    for phase=1:n_phases % loop over phase
        for reg=1:n_regions % loop over regions
            for t=1:n_types
                tmp.(types{t}) = [];
                %C_th.(types{t}) = [];
            end
            for subj=1:n_subjects % loop over subjects
                type = results_all{subj}.type;
                tmp.(type) = horzcat(tmp.(type),results_all{subj}.results.BLP.(phases{phase}).(regions{reg}).BLP); % concatenate data of the same type
                f = results_all{subj}.results.BLP.(phases{phase}).(regions{reg}).fc_envelope;
            end
    
            y = mean(tmp.(types{tp}),2); % take average over types
            y_err = std(tmp.(types{tp}),[],2)./sqrt(size(tmp.(types{tp}),2)); % get S.E over types
            x = f;
            BLP_stats_out.(phases{phase}).(regions{reg}).(types{tp}).coh_freq = f;
            BLP_stats_out.(phases{phase}).(regions{reg}).(types{tp}).coh_mean = y;
            BLP_stats_out.(phases{phase}).(regions{reg}).(types{tp}).coh_se = y_err;
            BLP_stats_out.(phases{phase}).(regions{reg}).(types{tp}).individual = tmp.(types{tp});
            subplot(1,n_regions,reg);xlabel('frequency (Hz)');title(regions_name{reg})
            plot_shadedErrorBar(x,y,y_err,'lineProps',{'color',phcolors.(phases{phase})});
            hold all
            axis square
            set_figure_properties();
        end
    end
    % create legend
    dim = [.9 .8 .2 .2];
    go_down = [0,0.035,0,0];
    for phase=1:n_phases
        str = phases{phase};
        dim_to_use = dim - (phase-1)*go_down;
        annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',phcolors.(phases{phase}), ...
            'FontSize',20);
    end
    subplot(1,n_regions,1);ylabel('BLP');
end
      
%% calculating normalized coherency and statistics
to_stats = struct;
for phase=2:n_phases % all phases except baseline
    for reg=1:n_regions % loop over regions
        for t=1:n_types
            tmp.(types{t}) = [];
        end
        for subj=1:n_subjects % for each subject
            type = results_all{subj}.type; % get the type
            tmp_coh = results_all{subj}.results.BLP.(phases{phase}).(regions{reg}).BLP;
            tmp_baseline = mean(results_all{subj}.results.BLP.(phases{1}).(regions{reg}).BLP,2); % get the average baseline
            tmp_subj_normalized = (tmp_coh - tmp_baseline)./((tmp_coh + tmp_baseline + eps)); % calculate modulation index
            tmp.(type) = cat(2,tmp.(type),tmp_subj_normalized); % concaenate modulation indices
            f = results_all{subj}.results.BLP.(phases{phase}).(regions{reg}).fc_envelope;
        end
        to_stats.(phases{phase}).(regions{reg}) = tmp;
        BLP_stats_out.(phases{phase}).(regions{reg}).modxfreq = tmp; %save modulation index per frequency
    end   
end
%Plot coherency modulation frequency resolved
BLP_ctrl = BLP_stats_out.(phases{2}).PFC_rs.modxfreq.ctrl;
BLP_exp = BLP_stats_out.(phases{2}).PFC_rs.modxfreq.exp;
BLP_ctrl_mean = mean(BLP_ctrl,2);
BLP_exp_mean = mean(BLP_exp,2);
BLP_ctrl_err = std(BLP_ctrl,[],2)/sqrt(size(BLP_ctrl,2));
BLP_exp_err = std(BLP_exp,[],2)/sqrt(size(BLP_exp,2));
figure('Name','Coherency Modulation  Frequency Resolved','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
plot_shadedErrorBar(f,BLP_ctrl_mean,BLP_ctrl_err,'lineProps',{'color',tycolors.ctrl}); 
hold on;
plot_shadedErrorBar(f,BLP_exp_mean,BLP_exp_err,'lineProps',{'color',tycolors.exp});
xlabel("frequency (Hz)");
ylabel("BLP modulation index");
dim = [.9 .8 .2 .2];
go_down = [0,0.035,0,0];
for ty = 1:length(types)
    str = types{ty};
    dim_to_use = dim - (ty-1)*go_down;
    annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',tycolors.(types{ty}), ...
             'FontSize',20);
end
set_figure_properties();

n_comp = (n_types*(n_types -1))/2; 
for phase=2:n_phases % loop over all phases except baseline
    figure('Name','BLP Modulation','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    St = sgtitle(phases{phase});
    St.FontSize =30;
    for reg=1:n_regions % loop over regions
        % init variables
        to_bar = zeros(n_types,1);
        err_bar = zeros(n_types,1);
        p_val = zeros(n_comp,1);
        to_test = cell(1,n_types);
        for type=1:n_types % loop over data types
            to_test{type} = to_stats.(phases{phase}).(regions{reg}).(types{type}); % get modulation indices of that phase,region,type
            to_bar(type,1) = mean(to_test{type}); % take the average
            err_bar(type,1) = std(to_test{type})/sqrt(length(to_test{type})); % take the S.E
        end
        comp = 1;
        for type1=1:n_types % loop over type1
            for type2=type1+1:n_types % loop over type2
                p_val(comp,1) = ranksum(to_test{type1},to_test{type2}); % compare the modulation indices between groups within region and phase
                comp = comp + 1;
            end
        end
        
        
        [~,~,~,p_val_corrected] = stat_fdr_bh(p_val); % correct p-value
        % plot results
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
        title(regions_name{reg});
        BLP_stats_out.compare.(phases{phase}).(regions{reg}).modulation_mean = to_bar;
        BLP_stats_out.compare.(phases{phase}).(regions{reg}).modulation_se = err_bar;
        BLP_stats_out.compare.(phases{phase}).(regions{reg}).modulation_p_val = p_val_corrected;
    end
    subplot(1,n_regions,1);ylabel('BLP modulation index');
    set_figure_properties()
end


%% comparing different phases
if results_all{1}.params.analyses.BLP.pooling == 0 % for no pooling do within group between phase analyses
    n_comp = (n_phases*(n_phases -1))/2; % get the number of comparisions
    for tp=1:n_types % loop over types
        figure('Name','Comparing BLP across conditions','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);sgtitle(types{tp});
        for reg=1:n_regions % loop over regions
            to_bar = zeros(n_phases,1);
            err_bar = zeros(n_phases,1);
            p_val = zeros(n_comp,1);
            to_test = cell(1,n_phases);
            for phase=1:n_phases % loop over phases
                for t=1:n_types % initialize dummy variable to concatenate data of the same type
                    tmp.(types{t}) = [];
                end
                for subj=1:n_subjects % loop over subjects
                    type = results_all{subj}.type; % get the subject type
                    tmp.(type) = horzcat(tmp.(type),results_all{subj}.results.BLP.(phases{phase}).(regions{reg}).BLP); % concatenate data of the same type
                end
                to_test{phase} = tmp.(types{tp}); % get the type you are analysing
                to_bar(phase,1) = mean(to_test{phase}); % get the average 
                err_bar(phase,1) = (std(to_test{phase})/sqrt(length(to_test{phase}))); % get S.E
            end
            comp = 1;
            for phase1=1:n_phases % loop over phase1 
                for phase2=phase1+1:n_phases % loop over phase2
                    [~,p_val(comp,1)] = ttest(to_test{phase1},to_test{phase2}); % pair test to see if the average has changed between phases
                    comp = comp + 1;
                end
            end
            
            [~,~,~,p_val_corrected] = stat_fdr_bh(p_val); % correct p-value
            % plot results
            subplot(1,n_regions,reg);
            bar(1:n_phases,to_bar);hold all;
            errorbar(1:n_phases,to_bar,err_bar,'.');
            ylim([0,0.5]);xticklabels(phases);
            comp = 1;
            for phase1=1:n_phases
                for phase2=phase1+1:n_phases
                    sigstar([phase1,phase2],p_val_corrected(comp));
                    comp = comp+1;
                end
            end
            title(regions_name{reg});
            BLP_stats_out.corr_vs.(regions{reg}).(types{tp}).mean = to_bar;
            BLP_stats_out.corr_vs.(regions{reg}).(types{tp}).se = err_bar;
            BLP_stats_out.corr_vs.(regions{reg}).(types{tp}).p_val = p_val_corrected;
            
            subplot(1,n_regions,1);ylabel('BLP correlation coeff');
            set_figure_properties()
        end
    end
else % if there was no pooling 
    to_stats = struct;
    for phase=2:n_phases % loop over all phases except baseline
        for reg=1:n_regions % loop over all regions
            for t=1:n_types % initialize variables for different groups (exp,ctrl,...)
                tmp.(types{t}) = [];
                tmp_BLP_all.(types{t}) = [];
                tmp_BLP_baseline.(types{t}) = [];
            end
            for subj=1:n_subjects % loop over subjects
                type = results_all{subj}.type; % get their type
                tmp_BLP = mean(results_all{subj}.results.BLP.(phases{phase}).(regions{reg}).BLP); % get the BLP of that phase
                tmp_baseline = results_all{subj}.results.BLP.(phases{1}).(regions{reg}).BLP; % get the average BLP of the baseline
                tmp_subj_normalized = (tmp_BLP - tmp_baseline)./((tmp_BLP + tmp_baseline + eps)); % calculate the modulation index
                tmp.(type) = cat(1,tmp.(type),tmp_subj_normalized'); % concatenate modulation indices of data of the same type
                tmp_BLP_all.(type) = cat(1,tmp_BLP_all.(type),tmp_BLP'); % ???
                tmp_BLP_baseline.(type) = cat(1,tmp_BLP_baseline.(type),tmp_baseline'); % ???
            end
            to_stats.(phases{phase}).(regions{reg}) = tmp;
            BLP_stats_out.compare.(phases{phase}).(regions{reg}).modulation_individual = tmp;
        end   
    end
    n_comp = (n_types*(n_types -1))/2; 
    for phase=2:n_phases % loop over all phases except baseline
        figure('Name','BLP Modulation','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
        St = sgtitle(phases{phase});
        St.FontSize =30;
        for reg=1:n_regions % loop over regions
            % init variables
            to_bar = zeros(n_types,1);
            err_bar = zeros(n_types,1);
            p_val = zeros(n_comp,1);
            to_test = cell(1,n_types);
            for type=1:n_types % loop over data types
                to_test{type} = to_stats.(phases{phase}).(regions{reg}).(types{type}); % get modulation indices of that phase,region,type
                to_bar(type,1) = mean(to_test{type}); % take the average
                err_bar(type,1) = std(to_test{type})/sqrt(length(to_test{type})); % take the S.E
            end
            comp = 1;
            for type1=1:n_types % loop over type1
                for type2=type1+1:n_types % loop over type2
                    p_val(comp,1) = ranksum(to_test{type1},to_test{type2}); % compare the modulation indices between groups within region and phase
                    comp = comp + 1;
                end
            end
            
            
            [~,~,~,p_val_corrected] = stat_fdr_bh(p_val); % correct p-value
            % plot results
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
            title(regions_name{reg});
            BLP_stats_out.compare.(phases{phase}).(regions{reg}).modulation_mean = to_bar;
            BLP_stats_out.compare.(phases{phase}).(regions{reg}).modulation_se = err_bar;
            BLP_stats_out.compare.(phases{phase}).(regions{reg}).modulation_p_val = p_val_corrected;
        end
        subplot(1,n_regions,1);ylabel('BLP modulation index');
        set_figure_properties()
    end
end


%% comparing different groups
n_comp = (n_types*(n_types -1))/2;
for phase=1:n_phases % loop over phases
    figure('Name','Comparing BLP across groups','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    st = sgtitle(phases{phase});
    st.FontSize =30;
    for reg=1:n_regions % loop over regions
        % init variables
        to_bar = zeros(n_types,1);
        err_bar = zeros(n_types,1);
        p_val = zeros(n_comp,1);
        to_test = cell(1,n_types);
        for t=1:n_types
                tmp.(types{t}) = [];
        end
        for subj=1:n_subjects % loop over subjects
            type = results_all{subj}.type; % get data type
            tmp.(type) = horzcat(tmp.(type),results_all{subj}.results.BLP.(phases{phase}).(regions{reg}).BLP); % concatenate data of same type
        end
        for type=1:n_types % loop over types
            to_test{type} = tmp.(types{type}); % get data of that type
            to_bar(type,1) = mean(to_test{type},'omitnan'); % get the average BLP for that type
            err_bar(type,1) = std(to_test{type},'omitnan')/sqrt(length(to_test{type})); % get the S.E 
        end
        comp = 1;
        for type1=1:n_types
            for type2=type1+1:n_types
                [~,p_val(comp,1)] = ttest2(to_test{type1},to_test{type2}); % compare data of the same phase and region across types
                comp = comp + 1; 
            end
        end
        
        
        [~,~,~,p_val_corrected] = stat_fdr_bh(p_val); % correct p-val
        % visualize
        subplot(1,n_regions,reg);
        bar(1:n_types,to_bar);hold all;
        errorbar(1:n_types,to_bar,err_bar,'.');
        ylim([0,0.5])
        xticklabels(types);
        set_figure_properties()
        comp = 1;
        for type1=1:n_types
            for type2=type1+1:n_types
                sigstar([type1,type2],p_val_corrected(comp));
                comp = comp + 1;
            end
        end
        title(regions_name{reg});
        BLP_stats_out.(phases{phase}).(regions{reg}).BLP_mean = to_bar;
        BLP_stats_out.(phases{phase}).(regions{reg}).BLP_se = err_bar;
        BLP_stats_out.(phases{phase}).(regions{reg}).BLP_p_val = p_val_corrected;
    end
    subplot(1,n_regions,1);ylabel('BLP correlation coeff');
    
end


end


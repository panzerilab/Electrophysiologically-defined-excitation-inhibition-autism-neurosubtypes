function PLV_stats_out = perform_PLV_stats(results_all,pl)
% Function to do PLV stats
% Inputs: 
% Results_all: a cell array of length of n_subjects iside each cell there 
% is a structure with the following fields:
% - type: the type of the subject (exp, ctrl,...)
% - params: a structure containing the parameters of all the analyses. for
% example params.PLV contains the parameters of the PLV analyses
% - results: a structure containing the results of various analyses. for
% example results.PLV contains the results of the PLV analyses
% pl: a structure containing the parameters for the group level 
% statistical analyses with the following fields:
% - sig_measure: a string to tell how you want to measure significancy of phase locking options are  {'PLV'} 
% or {'concentration'};
% - th a single value for the threshold for assessing the significance
% {value}>th ==> significant, for example, 0.1
% Outputs:
% PLV_stats_out: a structure containing various fields each with the
% results of the statistical analyses performed:
% - PLV_stats_out.(type_name).(phase_name).(region_str), is a structure
% containing various info for the for the data type (exp,ctrl,...) phase = phase_name
% (baseline,cno,..) and the region_str (PFC-rs,...)  with the following fields:
%   - phase_diff_pool: pooled distribution of all the angle differences. pooled over channels and subjects 
%   - p_val_uniform: the p-value of testing whether the distribution of angles 
%    in the phase_diff_pool is uniform or not
%   - mean_phase_different_from_zero: the p-value of testing whether the distribution of angles 
%    in the phase_diff_pool has a mean different than zero or not
%  - PLV_stats_out.phase_vs.(region_str).(type_name) contains following
%       fields:
%       - circ_mean: an array the length of n_phases which contains the
%       circular average angle of pooled data for specific data type in the
%       region_str for each phase in radians
%       - circ_se: an array the length of n_phases which contains the
%        standard error of the circular average angle of pooled data for specific data
%        type in the region_str for each phase in radians
%       - circ_pval: an array containing the p-values for statistical
%       comparision of the angle means of different phases within the same
%       group and region between different phases
%  - PLV_stats_out.compare.(phase_name).(region_str) contains following
%       fields for the phase = phase_name (except the baseline) and the
%       region = region_str
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
% -----------------------------------------------------------------------------
% this routine, loops over subjects, determines their type. concatenates
% the angle differences for the data of the same type. checks for the
% thresold, and for the data that exceeds the threshold (plv/concentaration is higher than
% certain value) plots the polar histogram of the angle differences for each phase and region.
% checks for the non-uniformity of the distribution using Rayleigh test (circ_rtest) and 
% checks whether the mean is different from zero or not using circ_mtest. 
% then within same data type, for each region, it checks whether the
% average angle has changed from one phase to another using Watson-Williams multi-sample test (circ_wwtest)
% followed by FDR correction and visualizes the results. After that for each data type, within each
% region and for all the phases except the baseline, it calculates a
% modulation index as follow: 
% MI = (PLV(region_xy,phase_z) - PLV(region_xy,baseline)/(PLV(region_xy,phase_z) + PLV(region_xy,baseline)
% and then compares the MIs across different data types using ranksum test
% followed by FDR correction.
% note: it will always consider the first phase as baseline
% note: LFP PLV is always measure between two regions, so by region we
% mean region-pair (for example PFC-rs)
% Requires circ_stat toolbox

%% initialize
n_subjects = length(results_all);
phases = fieldnames(results_all{1}.results.PLV);
regions = fieldnames(results_all{1}.results.PLV.(phases{1}));



type_all = cell(1,n_subjects);
for subj=1:n_subjects
    type_all{subj} = results_all{subj}.type;
end

types = unique(type_all,'stable');



if ~strcmp(types{1},'ctrl')
    types{2} = types{1};
    types{1} = 'ctrl'
end

phcolors.baseline = '#1AFF1A';
phcolors.cno=  '#FF00FF';

n_types = length(types);
n_subjects = length(results_all);
n_phases = length(phases);
n_regions = length(regions);

regions_name = regions;
for reg=1:n_regions
    temp_reg = regions{reg};
    idx_underline = strfind(temp_reg,'_');
    temp_reg(idx_underline) = '-';
    regions_name{reg} = temp_reg;
end
%% Plotting the histograms and showing which area is ahead of the other in phase
edges = 0:pi/16:2*pi;
for tp=1:n_types
    figure('Name','Histogram of phases','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    st = sgtitle(types{tp});
    st.FontSize = 20;
    count = 1;
    for phase=1:n_phases
        for reg=1:n_regions
            for t=1:n_types
                tmp.(types{t}) = [];
                tmp2.(types{t}) = [];
                tmp3.(types{t}) = [];
            end
            for subj=1:n_subjects
                type = results_all{subj}.type;
                tmp.(type) = horzcat(tmp.(type),results_all{subj}.results.PLV.(phases{phase}).(regions{reg}).phase_diff_mean);
                tmp2.(type) = horzcat(tmp2.(type),results_all{subj}.results.PLV.(phases{phase}).(regions{reg}).PLV);
                tmp3.(type) = horzcat(tmp3.(type),results_all{subj}.results.PLV.(phases{phase}).(regions{reg}).concentration);
            end
            to_hist = tmp.(types{tp});
            %Exclude phases whose PLV is less than 0.1
            if strcmp(pl.sig_measure,'PLV')
                to_hist_check = tmp2.(types{tp});
            elseif strcmp(pl.sig_measure,'concentration')
                to_hist_check = tmp3.(types{tp});
            else
                error('the requested significance assessment is not supported')
            end
            to_hist_check_th = pl.th;           
            to_hist = to_hist(to_hist_check>to_hist_check_th);
            if isempty(to_hist)
                to_hist = 0;
            end
            PLV_stats_out.(types{tp}).(phases{phase}).(regions{reg}).phase_diff_pool = wrapTo360(rad2deg(to_hist));
            PLV_stats_out.(types{tp}).(phases{phase}).(regions{reg}).p_val_uniform = circ_rtest(to_hist);
            PLV_stats_out.(types{tp}).(phases{phase}).(regions{reg}).mean_phase_different_from_zero = circ_mtest(to_hist,0);
            
            %subplot(n_phases,n_regions,count);
            count = count  +1;
            mean_angle = wrapTo2Pi(circ_mean(to_hist'));
            r = circ_r(to_hist');
            r = 0:0.01:r*0.5;
            
            str_title{1} = phases{phase};
            str_title{2} = regions_name{reg};
            str_title{3} = strcat('non-uniformity: ',getpstars(circ_rtest(to_hist)));
            if circ_mtest(to_hist,0)
                str_title{4} = ' significant different from zero';
            else
                str_title{4} = 'non-significant different zero';
            end
            polarhistogram(to_hist,edges,'normalization','probability','FaceColor',phcolors.(phases{phase}),'FaceAlpha',.3);
            %st = title(str_title);
            %st.FontSize = 20;
            %rlim([0,0.5]);
            hold all
            rl = rlim
            r = 0:0.01:rl(2);
            polarplot(mean_angle.*ones(size(r)),r,'linewidth',4,'color',phcolors.(phases{phase}));
            if phase ==1
                dim = [0.08,0.5,0.1,0.1];
            else
                dim = [0.8,0.5,0.1,0.1];
            end
            annotation('textbox',dim,'String',str_title,'FitBoxToText','on', 'FontSize',18);
        end

        dim = [.9 .8 .2 .2];
        go_down = [0,0.035,0,0];
        for phase=1:n_phases
            str = phases{phase};
            dim_to_use = dim - (phase-1)*go_down;
            annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',phcolors.(phases{phase}), ...
                'FontSize',20);
        end
    end   
end

%% comparing different phases
n_comp = (n_phases*(n_phases -1))/2;
for type=1:n_types
     figure('Name','Comparing phases','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
     st = sgtitle(types{type});
     st.FontSize = 30;
    for reg=1:n_regions
        to_bar = zeros(n_phases,1);
        err_bar = zeros(n_phases,1);
        p_val = zeros(n_comp,1);
        to_test = cell(1,n_phases);
        for phase=1:n_phases
            to_test{phase} = deg2rad(PLV_stats_out.(types{type}).(phases{phase}).(regions{reg}).phase_diff_pool');
            to_bar(phase,1) = rad2deg(circ_mean(to_test{phase}));
            err_bar(phase,1) = rad2deg(circ_confmean(to_test{phase}));
        end
        
        
        comp = 1;
        for phase1=1:n_phases
            for phase2=phase1+1:n_phases
                p_val(comp,1) = circ_wwtest(to_test{phase1},to_test{phase2});
                comp = comp + 1;
            end
        end
        
        [~,~,~,p_val_corrected] = stat_fdr_bh(p_val);
        
        subplot(1,n_regions,reg);
        bar(1:n_phases,to_bar);hold all;
        errorbar(1:n_phases,to_bar,err_bar,'.');
        ylim([-180,180]);xticklabels(phases);
        comp = 1;
        for phase1=1:n_phases
            for phase2=phase1+1:n_types
                sigstar([phase1,phase2],p_val_corrected(comp));
                comp = comp+1;
            end
        end
        title(regions_name{reg});
        PLV_stats_out.phase_vs.(regions{reg}).(types{type}).circ_mean = to_bar;
        PLV_stats_out.phase_vs.(regions{reg}).(types{type}).circ_se = err_bar;
        PLV_stats_out.phase_vs.(regions{reg}).(types{type}).circ_p_val = p_val_corrected;
        
        subplot(1,n_regions,1);ylabel('circular mean angle (degree)');
        set_figure_properties()
    end
end



%% calculating normalized PLV and statistics

to_stats = struct;
for phase=2:n_phases
    for reg=1:n_regions
        for t=1:n_types
            tmp.(types{t}) = [];
            tmp_PLV_all.(types{t}) = [];
            tmp_PLV_baseline.(types{t}) = [];
        end
        for subj=1:n_subjects
            type = results_all{subj}.type;
            tmp_PLV = results_all{subj}.results.PLV.(phases{phase}).(regions{reg}).PLV;
            tmp_baseline = results_all{subj}.results.PLV.(phases{1}).(regions{reg}).PLV;
            tmp_subj_normalized = (tmp_PLV - tmp_baseline)./((tmp_PLV + tmp_baseline + eps));
            tmp.(type) = cat(1,tmp.(type),tmp_subj_normalized');
            tmp_PLV_all.(type) = cat(1,tmp_PLV_all.(type),tmp_PLV');
            tmp_PLV_baseline.(type) = cat(1,tmp_PLV_baseline.(type),tmp_baseline');
        end
        to_stats.(phases{phase}).(regions{reg}) = tmp;
        PLV_stats_out.compare.(phases{phase}).(regions{reg}).modulation_individual = tmp;
        for tp=1:n_types 
            PLV_stats_out.(types{tp}).(phases{phase}).(regions{reg}).PLV_individual = tmp_PLV_all.(types{tp});
            PLV_stats_out.(types{tp}).baseline.(regions{reg}).PLV_individual = tmp_PLV_baseline.(types{tp});
        end
    end
    
end


n_comp = (n_types*(n_types -1))/2;
for phase=2:n_phases
    figure('Name','PLV Modulation','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    st = sgtitle(phases{phase});
    st.FontSize = 30;
    for reg=1:n_regions
        to_bar = zeros(n_types,1);
        err_bar = zeros(n_types,1);
        p_val = zeros(n_comp,1);
        to_test = cell(1,n_types);
        for type=1:n_types
            to_test{type} = to_stats.(phases{phase}).(regions{reg}).(types{type});
            to_bar(type,1) = mean(to_test{type});
            err_bar(type,1) = std(to_test{type})/sqrt(length(to_test{type}));
        end
        comp = 1;
        for type1=1:n_types
            for type2=type1+1:n_types
                p_val(comp,1) = ranksum(to_test{type1},to_test{type2});
                comp = comp + 1;
            end
        end
        
        
        [~,~,~,p_val_corrected] = stat_fdr_bh(p_val);
           
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
        PLV_stats_out.compare.(phases{phase}).(regions{reg}).modulation_mean = to_bar;
        PLV_stats_out.compare.(phases{phase}).(regions{reg}).modulation_se = err_bar;
        PLV_stats_out.compare.(phases{phase}).(regions{reg}).modulation_p_val = p_val_corrected;
    end
    subplot(1,n_regions,1);ylabel('PLV modulation index');
    set_figure_properties()
end

end





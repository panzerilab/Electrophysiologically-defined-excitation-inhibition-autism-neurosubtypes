function coherency_stats_out = perform_coherency_stats(results_all,coh)
% Function to do coherency stats
% Inputs: 
% Results_all: a cell array of length of n_subjects iside each cell there 
% is a structure with the following fields:
% - type: the type of the subject (exp, ctrl,...)
% - params: a structure containing the parameters of all the analyses. for
% example params.coherency contains the parameters of the coherency analyses
% - results: a structure containing the results of various analyses. for
% example results.coherency contains the results of the coherency analyses
% coh: a structure containing the parameters for the group level 
% statistical analyses with the following fields:
% - bands: a cell array of requested length within each cell an array of 
% lenght two of positive numbers ([f1,f2], f2>f1) to determine the bands
% for example {[0.1,1],[1,4],[4,8],[8,12],[12,30],[30,70]};
% - bands_name: a cell array with the same lenght of bands within each cell
% a string that determines the corresponding bands name for example: 
% {'infraslow','delta','theta','alpha','beta','gamma'};
% Outputs:
% coherency_stats_out: a structure containing various fields each with the
% results of the statistical analyses performed:
% - coherency_stats_out.(phase_name).(region_str).(type_name), is a structure
% containing various info for the data in phase = phase_name
% (baseline,cno,..) and the region_str (PFC-rs,...) for the data type
% (exp,ctrl,...) with the following fields:
%   - coh_freq: the frequency at which the coherency was calculated
%   - coh_mean: the average coherency
%   - coh_se: the standard error of the average of coherency
%   - individual: a matrix that at each column, the individual
%   coherency for an independent sample are stored. note that
%   independent sample does not necessary mean different subject. if there
%   was pooling or channel level analyses, each channel and each pooled
%   segment will be treated as independent
%  -  coherency_stats_out.(phase_name).(region_str) additionally contains following
%  fields for all the pahses except baseline
%       - modulation_individual: a structure wih fieldnames corresponding
%       to bands_name, inside each a cell of the length of n_types, within
%       each cell there are calculated individual modulation indices for
%       that data type. the types are sorted based on the visualization on
%       the figure, i.e. if exp is depicted first, the first element of
%       cell is for exp data
%       - modulation_mean: a matrix of size n_types x n_bands wich stores
%       the average modulation index of that region and phase. the order of
%       the rows are sorted based on the visualization on
%       the figure, i.e. if exp is depicted first, the first row of the 
%       matrix is for exp data
%       - modulation_se: a matrix of size n_types x n_bands wich stores
%       the standard error of the average modulation index of that region and 
%       phase. the order of the rows are sorted based on the visualization 
%       on the figures, i.e. if exp is depicted first, the first row of the 
%       matrix is for exp data
%       - modulation_p_val: a matrix of size n_comparasions x n_bands which
%       sores the p-values of comparing modulation indices between differnt
%       data types.
% -----------------------------------------------------------------------------
% this routine, loops over subjects, determines their type. concatenates
% the coherency of independent samples from subjects of the same type. then 
% for each phase,region visualizes the avg +/- se of the coherency 
% on a seperate figure for each data type. it proceeds to calculating the
% modulation index for each independent sample for each phase except the baseline.
% within each region as follows: 
% MI = (coherency(f,region_xy,phase_y) - coherency(f,region_xy,baseline)/(coherency(f,region_xy,phase_y) +coherency(f,region_xy,baseline)
% then for each band defined by the user. it takes the average modulation index
% within that frequency band for each independent sample and compares them
% for each region between different dataypes using ranksum test followed by
% FDR correction. 
% note: it will always consider the first phase as baseline
% note: coherency is always measure between two regions, so by region we
% mean region-pair (for example PFC-rs)
%% initialize
bands = coh.bands;
bands_name = coh.bands_name;
n_bands = length(bands);
flim = coh.flim;


n_subjects = length(results_all);
phases = fieldnames(results_all{1}.results.coherency);
regions = fieldnames(results_all{1}.results.coherency.(phases{1}));

n_phases = length(phases);
n_regions = length(regions);



regions_name = regions;
% convert region_str from reg1_reg2 to reg1-reg2
for reg=1:n_regions
    temp_reg = regions{reg};
    idx_underline = strfind(temp_reg,'_');
    temp_reg(idx_underline) = '-';
    regions_name{reg} = temp_reg;
end

type_all = cell(1,n_subjects);
for subj=1:n_subjects
    type_all{subj} = results_all{subj}.type;
end

types = unique(type_all,'stable'); % get how many types you have
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
if isfield(results_all{1}.results.coherency.(phases{1}).(regions{1}),'coh_th')
 shuffle = true
else
 shuffle = false
end
for tp=1:n_types % loop over types
    figure('Name','Raw Coherency','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    st = sgtitle(types{tp});
    st.FontSize = 30;
    for phase=1:n_phases % loop over phase
        for reg=1:n_regions % loop over regions
            for t=1:n_types
                tmp.(types{t}) = [];
                C_th.(types{t}) = [];
            end
            for subj=1:n_subjects % loop over subjects
                type = results_all{subj}.type;
                tmp.(type) = horzcat(tmp.(type),results_all{subj}.results.coherency.(phases{phase}).(regions{reg}).coherency); % concatenate data of the same type
                f = results_all{subj}.results.coherency.(phases{phase}).(regions{reg}).freq;
                if shuffle
                    C_th.(type) = horzcat(C_th.(type),results_all{subj}.results.coherency.(phases{phase}).(regions{reg}).coh_th) ;
                end
            end
    
            y = mean(tmp.(types{tp}),2); % take average over types
            y_err = std(tmp.(types{tp}),[],2)./sqrt(size(tmp.(types{tp}),2)); % get S.E over types
            x = f;
            coherency_stats_out.(phases{phase}).(regions{reg}).(types{tp}).coh_freq = f;
            coherency_stats_out.(phases{phase}).(regions{reg}).(types{tp}).coh_mean = y;
            coherency_stats_out.(phases{phase}).(regions{reg}).(types{tp}).coh_se = y_err;
            coherency_stats_out.(phases{phase}).(regions{reg}).(types{tp}).individual = tmp.(types{tp});
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
    if shuffle
        C_sl = mean(C_th.(types{tp}));
        for phase=1:n_phases 
            coherency_stats_out.(phases{phase}).C_sl = C_sl;
        end
    else
        %Compute significance level: analytical
        dof = 2*60*1000/length(f);
        q = dof/2 -1;
        C_sl = 1 - 0.05.^(1/q);
    end
    
    yline(C_sl,'-','0.05','LineWidth',2);
    subplot(1,n_regions,1);ylabel('coherency');
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
            tmp_coh = results_all{subj}.results.coherency.(phases{phase}).(regions{reg}).coherency;
            tmp_baseline = mean(results_all{subj}.results.coherency.(phases{1}).(regions{reg}).coherency,2); % get the average baseline
            tmp_subj_normalized = (tmp_coh - tmp_baseline)./((tmp_coh + tmp_baseline + eps)); % calculate modulation index
            tmp.(type) = cat(2,tmp.(type),tmp_subj_normalized); % concaenate modulation indices
            f = results_all{subj}.results.coherency.(phases{phase}).(regions{reg}).freq;
        end
        to_stats.(phases{phase}).(regions{reg}) = tmp;
        coherency_stats_out.(phases{phase}).(regions{reg}).modxfreq = tmp; %save modulation index per frequency
    end
   
end

n_comp = (n_types*(n_types -1))/2; % calculate the number of comparisions
for phase=2:n_phases % loop over phases except baseline
    figure('Name','Coherency Modulation','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    st = sgtitle(phases{phase});
    st.FontSize = 30;
    for reg=1:n_regions % loop over regions
        to_bar = zeros(n_types,n_bands);
        err_bar = zeros(n_types,n_bands);
        p_val = zeros(n_comp,n_bands);
        for band=1:n_bands % loop over bands
            idx_band = f >= bands{band}(1) & f < bands{band}(2); % get the indices of the band
            to_test = cell(1,n_types);
            for type=1:n_types % for each type
                tmp_band = to_stats.(phases{phase}).(regions{reg}).(types{type})(idx_band,:); % get the modulation index of the band
                to_test{type} = mean(tmp_band); % average over frequency
                to_bar(type,band) = mean(to_test{type}); % average over samples
                err_bar(type,band) = std(to_test{type})/sqrt(length(to_test{type})); % S.E over samples
            end
            comp = 1;
            for type1=1:n_types
                for type2=type1+1:n_types
                    p_val(comp,band) = ranksum(to_test{type1},to_test{type2}); % compare data type 1 and data type 2 
                    comp = comp + 1;
                end
            end         
            coherency_stats_out.(phases{phase}).(regions{reg}).modulation_individual.(bands_name{band}) = to_test;
        end
        [~,~,~,p_val_corrected] = stat_fdr_bh(p_val); % correct p-values
        subplot(1,n_regions,reg);barsem(to_bar,err_bar,p_val_corrected);ylim([-1,1]);
        xlabel('bands');xticklabels(bands_name);title(regions_name{reg});
        set_figure_properties();
        coherency_stats_out.(phases{phase}).(regions{reg}).modulation_mean = to_bar;
        coherency_stats_out.(phases{phase}).(regions{reg}).modulation_se = err_bar;
        coherency_stats_out.(phases{phase}).(regions{reg}).modulation_p_val = p_val_corrected;
    end
    % create a legend
    dim = [.9 .8 .2 .2];
    go_down = [0,0.035,0,0];
    for type=1:n_types
        str = types{type};
        dim_to_use = dim - (type-1)*go_down;
        annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',tycolors.(types{type}), ...
            'FontSize',20);
    end
    subplot(1,n_regions,1);ylabel('coherency modulation index');
end
%Plot coherency modulation frequency resolved and compute significative
%clusters
for phase=2:n_phases
    for reg=1:n_regions % loop over regions
        figure('Name','Coherency Modulation  (Frequency Resolved)','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
        st = sgtitle(phases{phase});
        st.FontSize = 30;
        to_test = cell(1,n_types);
        for type=1:n_types % for each type
            to_test{type} = coherency_stats_out.(phases{phase}).(regions{reg}).modxfreq.(types{type});
            x = f;
            y = mean(to_test{type},2);
            y_err = std(to_test{type},[],2)/sqrt(size(to_test{type},2));
            plot_shadedErrorBar(x,y,y_err,'lineProps',{'color',tycolors.(types{type})}); 
        end
        xlabel("frequency (Hz)");
        ylabel("coherency modulation index");
        dim = [.9 .8 .2 .2];
        go_down = [0,0.035,0,0];
        for ty = 1:length(types)
            str = types{ty};
            dim_to_use = dim - (ty-1)*go_down;
            annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',tycolors.(types{ty}), ...
                     'FontSize',20);
        end
        ilim = find(x>flim,1);
        if isempty(ilim)
            ilim = length(f)
        end
        set_figure_properties();
        %Perform cluster correction analyses and shade significative
        %regions in the plot
        hold on;
        [clusters,p_values] = cluster_correction(to_test{1},to_test{2},f,ilim);
        coherency_stats_out.(phases{phase}).(regions{reg}).clusters = clusters;
        coherency_stats_out.(phases{phase}).(regions{reg}).clusters_p_val = p_values;
    end
end


end


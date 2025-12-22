function spectrum_stats_out = perform_spectrum_stats(results_all,spec)
% Function to do spectrum stats
% Inputs: 
% Results_all: a cell array of length of n_subjects iside each cell there 
% is a structure with the following fields:
% - type: the type of the subject (exp, ctrl,...)
% - params: a structure containing the parameters of all the analyses. for
% example params.spectrum contains the parameters of the spectrum analyses
% - results: a structure containing the results of various analyses. for
% example results.coherency contains the results of the coherency analyses
% spec: a structure containing the parameters for the group level 
% statistical analyses with the following fields:
% - bands: a cell array of requested length within each cell an array of 
% lenght two of positive numbers ([f1,f2], f2>f1) to determine the bands
% for example {[0.1,1],[1,4],[4,8],[8,12],[12,30],[30,70]};
% - bands_name: a cell array with the same lenght of bands within each cell
% a string that determines the corresponding bands name for example: 
% {'infraslow','delta','theta','alpha','beta','gamma'};
% - fof: 0/1 if you don't want (want) the fof analyses
% outputs:
% spectrum_stats_out: a structure containing various fields each with the
% results of the statistical analyses performed:
% - spectrum_stats_out.(phase_name).(region_name).(type_name), is a structure
% containing various info for the data in phase = phase_name
% (baseline,cno,..) and the region_name (PFC,rs,...) for the data type
% (exp,ctrl,...) with the following fields:
%   - spec_freq: the frequency at which the spectrum was calculated
%   - spec_mean: the average spectrum in log10
%   - spec_se: the standard error of the average of spectrum in log10 scale
%   - spec_individual: a matrix that at each column, the individual
%   spectrum (in log10) for an independent sample are stored. note that
%   independent sample does not necessary mean different subject. if there
%   was pooling or channel level analyses, each channel and each pooled
%   segment will be treated as independent
%  -  spectrum_stats_out.(phase_name).(region_name) additionally contains following
%  fields for all the pahses except baseline
%       - mod_xband: a structure wih fieldnames corresponding
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
% - spectrum_stats_out.(phase_name).(region_name) will contain this
% additional fields if the fof analyses is requestd:
%  -fof_individual: a structure with fieldnames of different types (i.e.
%  exp, ctrl) inside each fof results are stored for each independent
%  sample for each band
%  -fof_mean: a structure with fieldnames of different types (i.e.
%  exp, ctrl) inside each the average fof results are stored for each band are stored    
%  -fof_err: a structure with fieldnames of different types (i.e.
%  exp, ctrl) inside each the standard error of average fof results are 
%  stored for each band are stored
% -----------------------------------------------------------------------------
% this routine, loops over subjects, determines their type. concatenates
% the spectrum of independent samples from subjects of the same type. then 
% for each phase,region visualizes the avg +/- se of the %spectrum in log10 
% on a seperate figure for each data type. it proceeds to calculating the
% modulation index for each independent sample for each phase except the baseline.
% within each region as follows: 
% MI = (power(f,region_x,phase_y) - power(f,region_x,baseline)/(power(f,region_x,phase_y) + power(f,region_x,baseline)
% then for each band defined by the user. it takes the average modulation index
% within that frequency band for each independent sample and cpmpares them
% for each region between different dataypes using ranksum test followed by
% FDR correction. 
% Additionally, if fof analyses is requested, it proceeds to calculate the
% fraction of power in each band over total power within each region and phase. visulizes it for each
% data type on a separete figure and compares their average across data types.
% if no pooling was requested, it will also do a cpmparision within data
% types across different phases using ranksum at each band followed by FDR
% note: it will always consider the first phase as baseline
%% initialzie variables

bands = spec.bands;
bands_name = spec.bands_name;
n_bands = length(bands);
flim = spec.flim;
n_subjects = length(results_all);
phases = fieldnames(results_all{1}.results.spectrum);
regions = fieldnames(results_all{1}.results.spectrum.(phases{1}));

n_phases = length(phases);
n_regions = length(regions);

type_all = cell(1,n_subjects);
for subj=1:n_subjects % get all the subjects types
    type_all{subj} = results_all{subj}.type;
end

types = unique(type_all,'stable'); % get how many types are there
if ~strcmp(types{1},'ctrl')
    types{2} = types{1};
    types{1} = 'ctrl';
end

n_types = length(types);

%% plotting spectrum
line_props = {'-b','-r','-g','-k'}; % colors for plotting 
tycolors.ctrl = "#0072BD";
tycolors.exp = "#D95319";
phcolors.baseline = '#1AFF1A';
phcolors.cno=  '#FF00FF';
phcolors.CNO=  '#FF00FF';

for tp=1:n_types % loop over groups
    figure('Name','Spectrum','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    st = sgtitle(types{tp});
    st.FontSize = 30;
    for phase=1:n_phases % loop over phases
        for reg=1:n_regions % loop over regions
            for t=1:n_types
                tmp.(types{t}) = [];
            end
            
            for subj=1:n_subjects % loop over subjects
                type = results_all{subj}.type; % get data type
                tmp.(type) = horzcat(tmp.(type),results_all{subj}.results.spectrum.(phases{phase}).(regions{reg}).spectrum); % concatanate data of the same type
                f = results_all{subj}.results.spectrum.(phases{phase}).(regions{reg}).spectrum_freq;
            end
            
            y = mean(log10(tmp.(types{tp})),2); % get average of the type in log10
            y_err = std(log10(tmp.(types{tp})),[],2)./sqrt(size(tmp.(types{tp}),2)); % get S.E of the type
            x = f;
            spectrum_stats_out.(phases{phase}).(regions{reg}).(types{tp}).spec_freq = f;
            spectrum_stats_out.(phases{phase}).(regions{reg}).(types{tp}).spec_mean = y;
            spectrum_stats_out.(phases{phase}).(regions{reg}).(types{tp}).spec_se = y_err;
            spectrum_stats_out.(phases{phase}).(regions{reg}).(types{tp}).spec_individual = log10(tmp.(types{tp}));
            subplot(1,n_regions,reg);xlabel('frequency (Hz)');title(regions{reg})
            plot_shadedErrorBar(x,y,y_err,'lineProps',{'color',phcolors.(phases{phase})}) % visualize the average +/- se
            hold all
            axis square
            set_figure_properties();
        end
    end
    dim = [.9 .8 .2 .2];
    go_down = [0,0.035,0,0];
    for phase=1:n_phases
        str = phases{phase};
        dim_to_use = dim - (phase-1)*go_down;
        annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',phcolors.(phases{phase}), ...
            'FontSize',20);
    end
    subplot(1,n_regions,1);ylabel('power');
end


%% calculating normalized spectrum and statistics
to_stats = struct;
for phase=2:n_phases % for phases except baseline
    for reg=1:n_regions % for all regions
        for t=1:n_types
            tmp.(types{t}) = [];
        end
        for subj=1:n_subjects % loop over subject
            type = results_all{subj}.type; % get data type
            tmp_spec = results_all{subj}.results.spectrum.(phases{phase}).(regions{reg}).spectrum; % get spectrum
            tmp_baseline = mean(results_all{subj}.results.spectrum.(phases{1}).(regions{reg}).spectrum,2); % get average baseline spectrum
            tmp_subj_normalized = (tmp_spec - tmp_baseline)./((tmp_spec + tmp_baseline + eps)); % calculate modulation index
            tmp.(type) = cat(2,tmp.(type),tmp_subj_normalized); % concatenate
            f = results_all{subj}.results.spectrum.(phases{phase}).(regions{reg}).spectrum_freq;
        end
        to_stats.(phases{phase}).(regions{reg}) = tmp;  % store for later: modulation index x frequency
        spectrum_stats_out.(phases{phase}).(regions{reg}).mod_xfr = tmp;
    end

    
end


n_comp = (n_types*(n_types -1))/2; % calculate how many comp are needed
for phase=2:n_phases % loop over phases except baseline
    figure('Name','Spectrum Modulation','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    st = sgtitle(phases{phase});
    st.FontSize = 30;
    for reg=1:n_regions % loop over each region
        to_bar = zeros(n_types,n_bands);
        err_bar = zeros(n_types,n_bands);
        p_val = zeros(n_comp,n_bands);
        for band=1:n_bands % for each band
            idx_band = f >= bands{band}(1) & f < bands{band}(2); % find the indices of that band
            to_test = cell(1,n_types);
            for type=1:n_types % for each data type
                tmp_band = to_stats.(phases{phase}).(regions{reg}).(types{type})(idx_band,:); % get modulation indices of the band
                to_test{type} = mean(tmp_band); % take the mean in frequency
                to_bar(type,band) = mean(to_test{type}); % take the mean across samples
                err_bar(type,band) = std(to_test{type})/sqrt(length(to_test{type})); % take the S.E across samples
            end
            comp = 1;
            for type1=1:n_types % loop over group 1
                for type2=type1+1:n_types % loop over group 2
                    p_val(comp,band) = ranksum(to_test{type1},to_test{type2}); % do the test
                    comp = comp + 1;
                end
            end
            spectrum_stats_out.(phases{phase}).(regions{reg}).mod_xband.(bands_name{band}) = to_test;
        end
        [~,~,~,p_val_corrected] = stat_fdr_bh(p_val); % correct for multiple comparision
        subplot(1,n_regions,reg);barsem(to_bar,err_bar,p_val_corrected);ylim([-1,1]);
        xlabel('bands');xticklabels(bands_name);set_figure_properties();title(regions{reg});
        spectrum_stats_out.(phases{phase}).(regions{reg}).modulation_mean = to_bar;
        spectrum_stats_out.(phases{phase}).(regions{reg}).modulation_se = err_bar;
        spectrum_stats_out.(phases{phase}).(regions{reg}).modulation_p_val = p_val_corrected;
    end
    % put the legends
    dim = [.9 .8 .2 .2];
    go_down = [0,0.035,0,0];
    for type=1:n_types
        str = types{type};
        dim_to_use = dim - (type-1)*go_down;
        annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',tycolors.(types{type}), ...
            'FontSize',20);
    end
    subplot(1,n_regions,1);ylabel('spectrum modulation index');
end

%Plot spectrum modulation frequency resolved and compute significative
%clusters
for phase=2:n_phases
    for reg=1:n_regions % loop over regions
        figure('Name',' Spectrum Modulation (Frequency Resolved)','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
        st = sgtitle(regions{reg});
        st.FontSize = 30;
        to_test = cell(1,n_types);
        for type=1:n_types % for each type
            to_test{type} = spectrum_stats_out.(phases{phase}).(regions{reg}).mod_xfr.(types{type}) ;
            x = f;
            y = mean(to_test{type},2);
            y_err = std(to_test{type},[],2)/sqrt(size(to_test{type},2));
            plot_shadedErrorBar(x,y,y_err,'lineProps',{'color',tycolors.(types{type})}); 
        end
        xlabel("frequency (Hz)");
        ylabel("spectrum modulation index");
        dim = [.9 .8 .2 .2];
        go_down = [0,0.035,0,0];
        for ty = 1:length(types)
            str = types{ty};
            dim_to_use = dim - (ty-1)*go_down;
            annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',tycolors.(types{ty}), ...
                     'FontSize',20);
        end
        ilim = find(f>flim,1);
        if isempty(ilim)
            ilim = length(f);
        end
        set_figure_properties();
        %Perform cluster correction analyses and shade significative
        %regions in the plot
        hold on;
        [clusters,p_values] = cluster_correction(to_test{1},to_test{2},f,ilim);
        spectrum_stats_out.(phases{phase}).(regions{reg}).clusters = clusters;
        spectrum_stats_out.(phases{phase}).(regions{reg}).clusters_p_val = p_values;
    end
end




%% fof analyses
if spec.fof == 1 % only if fof is requested
    %% calculate the fof and compare across types
    for reg = 1:n_regions
        tmp_fof_baseline.(regions{reg}) = [];
    end
    
    for phase=1:n_phases % loop over phase
        figure('Name','Fraction of frequency','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
        st = sgtitle(phases{phase});
        st.FontSize = 30;
        for region=1:n_regions % loop over region
            for tp=1:n_types % loop over data types
                fof_all.(types{tp}) = []; % init the structure to restore fof
                fof_all_xfr.(types{tp}) = [];
                tmp_mi.(types{tp}) = [];
            end
            
            for subj=1:n_subjects % loop over subjects
                tmp = results_all{subj}.results.spectrum;
                type = results_all{subj}.type; % get the type of subject
                spect = tmp.(phases{phase}).(regions{region}).spectrum; % get the spectrum
                f = tmp.(phases{phase}).(regions{region}).spectrum_freq;
                n_trials = size(spect,2); % get the number of samples
                fof_bands = zeros(n_trials,n_bands);
                fof_freq = zeros(n_trials,length(f)); %for frquency resolved analysis
                for tr=1:n_trials % loop over samples
                    s = spect(:,tr);
                    fof_freq(tr,:) =  100*s/(sum(s));
                    for band=1:n_bands % loop over band
                        idx_band = f >= bands{band}(1) & f < bands{band}(2); % find the indices of the band
                        fof_bands(tr,band) = (sum(s(idx_band))/sum(s))*100; % calculate percentage of power in the band
                    end
                end
                fof_all.(type) = vertcat(fof_all.(type),fof_bands); % concatenate fof of the same type
                fof_all_xfr.(type) = vertcat(fof_all_xfr.(type),fof_freq); % concatenate fof of the same type
                %calculate modulation index of fof
                if phase ==1
                    % Save baseline for all subjects
                    tmp_fof_baseline.(regions{region}) = vertcat(tmp_fof_baseline.(regions{region}),mean(fof_freq,1));
                end
                if phase ==2 
                    fof_b = tmp_fof_baseline.(regions{region});
                    mod_index = (fof_freq - fof_b(subj,:))./ (fof_freq + fof_b(subj,:));
                    tmp_mi.(type) = vertcat(tmp_mi.(type),mod_index);
                end
            end
            spectrum_stats_out.(phases{phase}).(regions{region}).fof_individual_xband = fof_all;
            spectrum_stats_out.(phases{phase}).(regions{region}).fof_individual_xfr = fof_all_xfr;
            if phase ==2
                spectrum_stats_out.(phases{phase}).(regions{region}).fof_mod_xfr = tmp_mi;
            end
            to_bar = zeros(n_types,n_bands);
            err_bar = zeros(n_types,n_bands);
            for type=1:n_types % loop over each type
                
                to_bar(type,:) = mean(fof_all.(types{type})); % get the mean fof for the type
                spectrum_stats_out.(phases{phase}).(regions{region}).fof_mean.(types{type}) = mean(fof_all.(types{type}));
                err_bar(type,:) = std(fof_all.(types{type}))/sqrt(size(fof_all.(types{type}),1));
                spectrum_stats_out.(phases{phase}).(regions{region}).fof_err.(types{type}) = std(fof_all.(types{type}))/sqrt(size(fof_all.(types{type}),1));
            end
            
            n_comp = (n_types*(n_types -1))/2; % calculate the number of comps
            p_val = zeros(n_comp,n_bands);
            comp = 1;
            for type1=1:n_types
                for type2=type1+1:n_types
                    for band=1:n_bands % for each band
                        data1 = fof_all.(types{type1});
                        data2 = fof_all.(types{type2});
                        [~,p_val(comp,band)] = ttest2(data1(:,band),data2(:,band)); % compare type1 and type 2
                        
                    end
                    comp = comp + 1;
                end
            end
            [~,~,~,p_val_corrected] = stat_fdr_bh(p_val); % correct the p-vals
            subplot(1,n_regions,region);barsem(to_bar,err_bar,p_val_corrected);
            title(regions{region});xticklabels(bands_name)
            set_figure_properties()
        end
        % make the legend
        dim = [.9 .8 .2 .2];
        go_down = [0,0.035,0,0];
        for type=1:n_types
            str = types{type};
            dim_to_use = dim - (type-1)*go_down;
            annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',tycolors.(types{type}),'FontSize',20);
        end
        subplot(1,n_regions,1);ylabel('fof (%)');
    end
    % Plot FOF mod x frequency
    for phase=2:n_phases
        for reg=1:n_regions % loop over regions
            figure('Name','FoF  Modulation (Frequency Resolved)','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
            st = sgtitle(regions{reg});
            st.FontSize = 30;
            to_test = cell(1,n_types);
            for type=1:n_types % for each type
                to_test{type} = spectrum_stats_out.(phases{2}).(regions{reg}).fof_mod_xfr.(types{type})';
                x = f;
                y = mean(to_test{type},2);
                y_err = std(to_test{type},[],2)/sqrt(size(to_test{type},2));
                plot_shadedErrorBar(x,y,y_err,'lineProps',{'color',tycolors.(types{type})}); 
            end
            xlabel("frequency (Hz)");
            ylabel("fof modulation index");
            dim = [.9 .8 .2 .2];
            go_down = [0,0.035,0,0];
            for ty = 1:length(types)
                str = types{ty};
                dim_to_use = dim - (ty-1)*go_down;
                annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',tycolors.(types{ty}), ...
                         'FontSize',20);
            end
            ilim = find(f>flim,1);
            if isempty(ilim)
                ilim = length(f);
            end
            set_figure_properties();
            %Perform cluster correction analyses and shade significative
            %regions in the plot
            hold on;
            [clusters,p_values] = cluster_correction(to_test{1},to_test{2},f,ilim);
            spectrum_stats_out.(phases{phase}).(regions{reg}).clusters = clusters;
            spectrum_stats_out.(phases{phase}).(regions{reg}).clusters_p_val = p_values;
        end
    end
    
    %FOF mod x bands
    n_comp = (n_types*(n_types -1))/2; % calculate how many comp are needed
    for phase=2:n_phases % loop over phases except baseline
        figure('Name','FOF Modulation','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
        st = sgtitle(phases{phase});
        st.FontSize = 30;
        for reg=1:n_regions % loop over each region
            to_bar = zeros(n_types,n_bands);
            err_bar = zeros(n_types,n_bands);
            p_val = zeros(n_comp,n_bands);
            for band=1:n_bands % for each band
                idx_band = f >= bands{band}(1) & f < bands{band}(2); % find the indices of that band
                to_test = cell(1,n_types);
                for type=1:n_types % for each data type
                    tmp_band = spectrum_stats_out.(phases{phase}).(regions{reg}).fof_mod_xfr.(types{type})(:,idx_band); % get modulation indices of the band
                    to_test{type} = mean(tmp_band); % take the mean in frequency
                    to_bar(type,band) = mean(to_test{type}); % take the mean across samples
                    err_bar(type,band) = std(to_test{type})/sqrt(length(to_test{type})); % take the S.E across samples
                end
                comp = 1;
                for type1=1:n_types % loop over group 1
                    for type2=type1+1:n_types % loop over group 2
                        p_val(comp,band) = ranksum(to_test{type1},to_test{type2}); % do the test
                        comp = comp + 1;
                    end
                end
                spectrum_stats_out.(phases{phase}).(regions{reg}).fof_mod_xband.(bands_name{band}) = to_test;
            end
            [~,~,~,p_val_corrected] = stat_fdr_bh(p_val); % correct for multiple comparision
            subplot(1,n_regions,reg);barsem(to_bar,err_bar,p_val_corrected);ylim([-1,1]);
            xlabel('bands');xticklabels(bands_name);set_figure_properties();title(regions{reg});
        end
        % put the legends
        dim = [.9 .8 .2 .2];
        go_down = [0,0.035,0,0];
        for type=1:n_types
            str = types{type};
            dim_to_use = dim - (type-1)*go_down;
            annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',tycolors.(types{type}), ...
                'FontSize',20);
        end
        subplot(1,n_regions,1);ylabel('fof modulation index');
    end

    %% compare across phases within groups
    if results_all{1}.params.analyses.spectrum.pooling == 0 % only do within groups if there is no pooling
        n_comp = (n_phases*(n_phases -1))/2; 
        for tp=1:n_types
            figure('Name','Comparing fof across conditions','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
            st = sgtitle(types{tp});
            st.FontSize = 30;
            for reg=1:n_regions
                
                p_val = zeros(n_comp,n_bands);
                to_bar = zeros(n_phases,n_bands);
                err_bar = zeros(n_phases,n_bands);
                for phase=1:n_phases
                    to_bar(phase,:) = mean(spectrum_stats_out.(phases{phase}).(regions{reg}).fof_individual_xband.(types{tp}));
                    err_bar(phase,:) = std(spectrum_stats_out.(phases{phase}).(regions{reg}).fof_individual_xband.(types{tp}))/sqrt(size(spectrum_stats_out.(phases{phase}).(regions{reg}).fof_individual_xband.(types{tp}),1));
                end
                comp = 1;
                for phase1=1:n_phases
                    for phase2=phase1+1:n_phases
                        to_test{phase1} = spectrum_stats_out.(phases{phase1}).(regions{reg}).fof_individual_xband.(types{tp});
                        to_test{phase2} = spectrum_stats_out.(phases{phase2}).(regions{reg}).fof_individual_xband.(types{tp});
                        for band=1:n_bands
                            [~,p_val(comp,band)] = ttest(to_test{phase1}(:,band),to_test{phase2}(:,band));
                        end
                        comp = comp + 1;
                    end
                end
                
                [~,~,~,p_val_corrected] = stat_fdr_bh(p_val);
                
                subplot(1,n_regions,reg);
                subplot(1,n_regions,reg);barsem(to_bar,err_bar,p_val_corrected);
                title(regions{reg});xticklabels(bands_name)
                set_figure_properties()
                title(regions{reg});xticklabels(bands_name)
                subplot(1,n_regions,1);ylabel('fof (%)');
                set_figure_properties()
            end
            for phase=1:n_phases
                str = phases{phase};
                dim_to_use = dim - (phase-1)*go_down;
                annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',phcolors.(phases{phase}),'FontSize',22);
            end
        end
    end

end



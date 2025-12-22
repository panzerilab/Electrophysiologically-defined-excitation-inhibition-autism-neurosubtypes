function GC_stats_out = perform_GC_stats(results_all,cause)

bands = cause.bands;
bands_name = cause.bands_name;
n_bands = length(bands);

n_subjects = length(results_all);
phases = fieldnames(results_all{1}.results.GC);
regions_tmp = fieldnames(results_all{1}.results.GC.(phases{1}));
regions = regions_tmp(~contains(regions_tmp,'Rsquared'));


n_phases = length(phases);
n_regions = length(regions);

regions_name = regions;
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

types = unique(type_all,'stable');
n_types = length(types);
if ~strcmp(types{1},'ctrl')
    types{2} = types{1};
    types{1} = 'ctrl';
end


%% plotting raw GC
line_props = {'-b','-r','-g','-k'}; % colors for plotting 
tycolors.ctrl = "#0072BD";
tycolors.exp = "#D95319";
phcolors.baseline = '#1AFF1A';
phcolors.cno=  '#FF00FF';
phcolors.CNO=  '#FF00FF';
regstyle.(regions{1}) = '-';
regstyle.(regions{2}) = '-.';

figure('Name','Raw Granger Causality','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
st = sgtitle("Raw GC \color{blue}{ctrl} vs \color{red}{exp}");
st.FontSize = 30;
for tp=1:n_types
    for phase=1:n_phases
        for reg=1:n_regions           
            for t=1:n_types
                tmp.(types{t}) = [];
            end
            GC_sh = [];
            for subj=1:n_subjects
                type = results_all{subj}.type;
                tmp.(type) = horzcat(tmp.(type),results_all{subj}.results.GC.(phases{phase}).(regions{reg}).infreq);
                GC_sh = horzcat(GC_sh,results_all{subj}.results.GC.(phases{phase}).(regions{reg}).infreq_sh);
                f = results_all{subj}.results.GC.(phases{phase}).(regions{reg}).freq;
            end
            
            y = mean(tmp.(types{tp}),2);
            y_err = std(tmp.(types{tp}),[],2)./sqrt(size(tmp.(types{tp}),2));
            x = f;
            GC_stats_out.(phases{phase}).(regions{reg}).(types{tp}).GC_freq = f;
            GC_stats_out.(phases{phase}).(regions{reg}).(types{tp}).GC_mean = y;
            GC_stats_out.(phases{phase}).(regions{reg}).(types{tp}).GC_se = y_err;
            GC_stats_out.(phases{phase}).(regions{reg}).(types{tp}).GC_individual = tmp.(types{tp});
            %subplot(2,2*round(n_regions/2),(phase-1)*2 + tp);
            subplot(1,n_types,tp);
            xlabel('frequency (Hz)');
            %title(regions_name{reg},'Color',tycolors.(types{tp}));  
            title(types{tp},'Color',tycolors.(types{tp}));  
            plot_shadedErrorBar(x,y,y_err,'lineProps',{'color',phcolors.(phases{phase}),'LineStyle',regstyle.(regions{reg})});
            hold all
            plot(f,mean(GC_sh,2),'--m');
            axis square
            set_figure_properties();
        end
    end
    dim = [.9 .8 .2 .2];
    go_down = [0,0.035,0,0];
    for phase=1:n_phases
        str = phases{phase};
        dim_to_use = dim - (phase-1)*go_down;
        annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',phcolors.(phases{phase}),'FontSize',18);
    end
    for reg=1:n_regions
        str = [regstyle.(regions{reg}), ' ',regions{reg}]
        dim_to_use = dim - (1+reg)*go_down;
        annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','FontSize',18,'interpreter','none');
    end
    subplot(1,n_types,tp);xlabel('frequency (Hz)');ylabel('GC');
    %subplot(2,2*round(n_regions/2),(reg-1)*2 + tp);;ylabel('GC');
end


%% Computing GC modulation index
to_stats = struct;
for phase=2:n_phases
    for reg=1:n_regions
        for t=1:n_types
            tmp.(types{t}) = [];
        end
        for subj=1:n_subjects
            type = results_all{subj}.type;
            tmp_cause = results_all{subj}.results.GC.(phases{phase}).(regions{reg}).infreq;
            tmp_baseline = results_all{subj}.results.GC.(phases{1}).(regions{reg}).infreq;
            tmp_subj_normalized = (tmp_cause - tmp_baseline)./((tmp_cause + tmp_baseline + eps));
            tmp.(type) = cat(2,tmp.(type),tmp_subj_normalized);
            f = results_all{subj}.results.GC.(phases{phase}).(regions{reg}).freq;
        end
        to_stats.(phases{phase}).(regions{reg}) = tmp;
    end
   
end

n_comp = (n_types*(n_types -1))/2;
for phase=2:n_phases
    figure('Name','Granger Causality Modulation Index','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    st = sgtitle(phases{phase});
    st.FontSize = 30;
    for reg=1:n_regions
        to_bar = zeros(n_types,n_bands);
        err_bar = zeros(n_types,n_bands);
        p_val = zeros(n_comp,n_bands);
        for band=1:n_bands
            idx_band = f >= bands{band}(1) & f < bands{band}(2);
            to_test = cell(1,n_types);
            for type=1:n_types
                tmp_band = to_stats.(phases{phase}).(regions{reg}).(types{type})(idx_band,:);
                to_test{type} = mean(tmp_band);
                to_bar(type,band) = mean(to_test{type});
                err_bar(type,band) = std(to_test{type})/sqrt(length(to_test{type}));
            end
            comp = 1;
            for type1=1:n_types
                for type2=type1+1:n_types
                    p_val(comp,band) = ranksum(to_test{type1},to_test{type2});
                    comp = comp + 1;
                end
            end         
            GC_stats_out.(phases{phase}).(regions{reg}).modulation_individual.(bands_name{band}) = to_test;
        end
        [~,~,~,p_val_corrected] = stat_fdr_bh(p_val);
        subplot(2,round(n_regions/2),reg);barsem(to_bar,err_bar,p_val_corrected);ylim([-1,1]);
        xlabel('bands');xticklabels(bands_name);title(regions_name{reg});
        set_figure_properties()
        GC_stats_out.(phases{phase}).(regions{reg}).modulation_mean = to_bar;
        GC_stats_out.(phases{phase}).(regions{reg}).modulation_se = err_bar;
        GC_stats_out.(phases{phase}).(regions{reg}).modulation_p_val = p_val_corrected;
    end
    dim = [.9 .8 .2 .2];
    go_down = [0,0.035,0,0];
    for type=1:n_types
        str = types{type};
        dim_to_use = dim - (type-1)*go_down;
        annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',tycolors.(types{type}), ...
            'FontSize',22);
    end
    subplot(2,round(n_regions/2),1);ylabel('GC modulation index');
    subplot(2,round(n_regions/2),round(n_regions/2)+1);ylabel('GC modulation index');
end

for phase=2:n_phases
    figure('Name','Granger Causality Modulation Index (freq. resolved)','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    st = sgtitle(phases{phase});
    st.FontSize = 30;
    for reg=1:n_regions
        gc_ctrl = to_stats.(phases{phase}).(regions{reg}).ctrl;
        gc_exp = to_stats.(phases{phase}).(regions{reg}).exp;
        gc_ctrl_mean = mean(gc_ctrl,2);
        gc_exp_mean = mean(gc_exp,2);
        gc_ctrl_err = std(gc_ctrl,[],2)/sqrt(size(gc_ctrl,2));
        gc_exp_err = std(gc_exp,[],2)/sqrt(size(gc_exp,2));        
        subplot(2,round(n_regions/2),reg);
        plot_shadedErrorBar(f,gc_ctrl_mean,gc_ctrl_err,'lineProps',{'color',tycolors.ctrl}); 
        hold on;
        plot_shadedErrorBar(f,gc_exp_mean,gc_exp_err,'lineProps',{'color',tycolors.exp});       
        

        xlabel('f(Hz)');title(regions_name{reg});
        set_figure_properties()
        ylabel('GC modulation index');
    end
     dim = [.9 .8 .2 .2];
    go_down = [0,0.035,0,0];
    for type=1:n_types
        str = types{type};
        dim_to_use = dim - (type-1)*go_down;
        annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',tycolors.(types{type}), ...
            'FontSize',22);
    end
end
end

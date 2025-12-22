function coherence_vstime(stats_coherence,n_pools)
    line_props = {'-b','-r','-g','-k'};
    colors = {'blue','red','green','black'};
    phases =  fieldnames(stats_coherence);
    modulation_ind = stats_coherence.(phases{2}).PFC_rs.modulation_individual;
    bands = fieldnames(modulation_ind);
    figure('Name','Coherency Modulation in Time','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    for band = 1:length(bands)
        mod = modulation_ind.(bands{band});
        mod_1 = reshape(mod{1},n_pools,[]);
        mod_2 = reshape(mod{2},n_pools,[]);
        mod_1_err = std(mod_1,[],2)./sqrt(size(mod_1,2));
        mod_2_err = std(mod_2,[],2)./sqrt(size(mod_2,2));
        y1 = mean(mod_1,2);
        y2 = mean(mod_2,2);
        x = 1:n_pools;
        % create a legend
        dim = [.9 .8 .2 .2];
        go_down = [0,0.035,0,0];
        subplot(3,2,band);ylabel('coh. mod. index');
        %for type=1:n_types
        %    str = types{type};
        %    dim_to_use = dim - (type-1)*go_down;
        %    annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',colors{type},'FontSize',12);
        %end
        plot_shadedErrorBar(x,y1,mod_1_err,'lineProps',line_props(1));
        hold on
        plot_shadedErrorBar(x,y2,mod_2_err,'lineProps',line_props(2));
        xlabel('time(min)');title(bands{band});
        set_figure_properties();
    end
    dim = [.9 .8 .2 .2];
    go_down = [0,0.035,0,0];
    types = fieldnames(stats_coherence.baseline.PFC_rs)
    for type=1:length(types)
        str = types{type};
        dim_to_use = dim - (type-1)*go_down;
        annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',colors{type},'FontSize',22);
    end
end
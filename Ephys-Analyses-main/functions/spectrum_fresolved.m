function spectrum_fresolved(stats_spectrum)
    phases = fieldnames(stats_spectrum);
    types = {"ctrl","exp"};
    tycolors.ctrl = "#0072BD";
    tycolors.exp = "#D95319";
    f = stats_spectrum.(phases{1}).PFC.ctrl.spec_freq;
    spmod_ctrl = stats_spectrum.(phases{2}).PFC.modxfr.ctrl;
    spmod_exp = stats_spectrum.(phases{2}).PFC.modxfr.exp;
    spmod_ctrl_mean = mean(spmod_ctrl,2);
    spmod_exp_mean = mean(spmod_exp,2);
    spmod_ctrl_err = std(spmod_ctrl,[],2)/sqrt(size(spmod_ctrl,2));
    spmod_exp_err = std(spmod_exp,[],2)/sqrt(size(spmod_exp,2));
    %Plot results
    figure('Name','Coherency Modulation  Frequency Resolved','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    plot_shadedErrorBar(f,spmod_ctrl_mean,spmod_ctrl_err,'lineProps',{'color',tycolors.ctrl}); 
    hold on;
    plot_shadedErrorBar(f,spmod_exp_mean,spmod_exp_err,'lineProps',{'color',tycolors.exp});
    xlabel("frequency (Hz)");
    ylabel("spectrum modulation index");
    dim = [.9 .8 .2 .2];
    go_down = [0,0.035,0,0];
    for ty = 1:length(types)
        str = types{ty};
        dim_to_use = dim - (ty-1)*go_down;
        annotation('textbox',dim_to_use,'String',str,'FitBoxToText','on','Color',tycolors.(types{ty}), ...
                 'FontSize',22);
    end
    set_figure_properties();
end

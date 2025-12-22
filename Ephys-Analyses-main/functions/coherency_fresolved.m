function [coh_ctrl,coh_exp] = coherency_fresolved(stats_coherency)
    phases = fieldnames(stats_coherency);
    types = {"ctrl","exp"};
    tycolors.ctrl = "#0072BD";
    tycolors.exp = "#D95319";
    f = stats_coherency.(phases{1}).PFC_rs.ctrl.coh_freq;
    coh_ctrl = stats_coherency.(phases{2}).PFC_rs.modxfreq.ctrl;
    coh_exp = stats_coherency.(phases{2}).PFC_rs.modxfreq.exp;
    coh_ctrl_mean = mean(coh_ctrl,2);
    coh_exp_mean = mean(coh_exp,2);
    coh_ctrl_err = std(coh_ctrl,[],2)/sqrt(size(coh_ctrl,2));
    coh_exp_err = std(coh_exp,[],2)/sqrt(size(coh_exp,2));
    figure('Name','Coherency Modulation  Frequency Resolved','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    plot_shadedErrorBar(f,coh_ctrl_mean,coh_ctrl_err,'lineProps',{'color',tycolors.ctrl}); 
    hold on;
    plot_shadedErrorBar(f,coh_exp_mean,coh_exp_err,'lineProps',{'color',tycolors.exp});
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
    set_figure_properties();
end


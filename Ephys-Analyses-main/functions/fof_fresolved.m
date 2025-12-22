function [fof_mod_ctrl,fof_mod_exp] = fof_fresolved(stats_spectrum)
    phases = fieldnames(stats_spectrum);
    types = {"ctrl","exp"};
    tycolors.ctrl = "#0072BD";
    tycolors.exp = "#D95319";
    %Compute modulation index x fraciton of frequency
    f = stats_spectrum.(phases{1}).PFC.ctrl.spec_freq;
    fof_baseline_ctrl = mean(stats_spectrum.(phases{1}).PFC.fof_individual_xfr.ctrl,1);
    fof_baseline_exp = mean(stats_spectrum.(phases{1}).PFC.fof_individual_xfr.exp,1);
    fof_CNO_ctrl = stats_spectrum.(phases{2}).PFC.fof_individual_xfr.ctrl;
    fof_CNO_exp = stats_spectrum.(phases{2}).PFC.fof_individual_xfr.exp;
    fof_mod_ctrl = (fof_CNO_ctrl - fof_baseline_ctrl)./(fof_CNO_ctrl + fof_baseline_ctrl);
    fof_mod_exp = (fof_CNO_exp - fof_baseline_exp)./(fof_CNO_exp + fof_baseline_exp);

    fofmod_ctrl_mean = mean(fof_mod_ctrl,1);
    fofmod_exp_mean = mean(fof_mod_exp,1);
    fofmod_ctrl_err = std(fof_mod_ctrl,[],1)/sqrt(size(fof_mod_ctrl,1));
    fofmod_exp_err = std(fof_mod_exp,[],1)/sqrt(size(fof_mod_exp,1));
    %Plot results
    figure('Name','FoF Modulation  Frequency Resolved','NumberTitle','off','units','normalized','outerposition',[0 0 0.9 0.9]);
    plot_shadedErrorBar(f,fofmod_ctrl_mean,fofmod_ctrl_err,'lineProps',{'color',tycolors.ctrl}); 
    hold on;
    plot_shadedErrorBar(f,fofmod_exp_mean,fofmod_exp_err,'lineProps',{'color',tycolors.exp});
    xlabel("frequency (Hz)");
    ylabel("fof modulation index");
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

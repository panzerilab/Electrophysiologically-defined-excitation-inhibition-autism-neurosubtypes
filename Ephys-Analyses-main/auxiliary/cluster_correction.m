function [clusters, p_values] = cluster_correction(mod_ctrl,mod_exp,f,flim)
    %Perform cluster-based permutation analyses and shade significative
    %clusters in the current plot
    [clusters, p_values, t_sums, permutation_distribution ] = permutest(mod_ctrl(1:flim,:),mod_exp(1:flim,:), ...
            false,[],10^4,true);
    n_clusters = sum(p_values < 0.05);
    lims = ylim;
    for i =1:n_clusters
        H = area([f(clusters{i}(1)) f(clusters{i}(end))],[lims(1) lims(1)],lims(2));
        set(H(1),'FaceColor','k');
        alpha(.3);
        p_string = getpstars(p_values(i));
        posfig = gca().Position;
        dim_to_use = [posfig(1)+ posfig(3)*mean(f(clusters{i}))/f(end),posfig(2)+posfig(4)+0.01,0.05,0.05];
        anText = annotation('textbox',dim_to_use,'String',p_string,'FitBoxToText','on','FontSize',15);
    end
end



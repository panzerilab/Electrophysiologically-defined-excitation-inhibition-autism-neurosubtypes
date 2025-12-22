function spaghetti_plot(x,y,g,cm)


cmp = cool(length(unique(cm)));
groups = unique(cm,'stable');
colors = zeros(length(g),3);
for i=1:length(colors)
    if ~isnumeric(groups)
        colors(i,:) = cmp(strcmp(groups,cm{i}),:);
    else
        if ~isnan(cm(i))
            colors(i,:) = cmp(find(groups == cm(i),1,'first'),:);
        end
    end
end


[~,idx_group] = unique(g,'stable');
idx_group = [1;idx_group];


figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(idx_group)-1
    line(x(idx_group(i):idx_group(i+1)),y(idx_group(i):idx_group(i+1)),'color',colors(idx_group(i),:)');
    hold all
end

n_groups = length(groups);
color_unique = unique(colors,'rows','stable');
dim = [.9 .8 .2 .2];
go_down = [0,0.035,0,0];
if ~isnumeric(groups)
    for grp=1:n_groups
        dim_to_use = dim - (grp-1)*go_down;
        a = annotation('textbox',dim_to_use,'str',groups{grp},'color',color_unique(grp,:),'FitBoxToText','on');
        a.FontSize = 12;
    end
else
    colormap cool
    colorbar('Ticks',[0,1],'TickLabels',[min(cm),max(cm)])
end
end


function  set_figure_properties()

ax = gca;
ax.FontSize = 22;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2.0;
end

end


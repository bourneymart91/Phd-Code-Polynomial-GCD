if(SETTINGS.PLOT_GRAPHS)
    figure_name = sprintf('%s : Condition APF', mfilename);
    figure('name',figure_name)
    hold on
    plot(1:1:length(condition),(condition));
    title(figure_name)
    xlabel('iteration')
    ylabel('condition')
end

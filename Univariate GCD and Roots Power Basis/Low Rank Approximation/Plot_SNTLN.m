if(SETTINGS.PLOT_GRAPHS)
    
    figure_name = sprintf('%s - Residuals',mfilename);
    figure('name',figure_name)
    hold on
    plot(log10(condition),'-s');
    hold off
    
    
end
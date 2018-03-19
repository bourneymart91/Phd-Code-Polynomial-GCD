
if(SETTINGS.PLOT_GRAPHS)
    
    
    title = sprintf('%s - Residuals',mfilename());
    xlabel('Iteration Number')
    ylabel('log_{10} Residual')
    figure('name',title)
    hold on
    plot(log10(condition),'-s')
    hold off
    
end


if SETTINGS.PLOT_GRAPHS_LRA
    
    fig_name = sprintf('%s : condition',mfilename);
    figure('name',fig_name)
    hold on
    plot(log10(vCondition),'-s','DisplayName','Condition')
    xlabel('iteration');
    ylabel('log_{10} condition');
    hold off
    
    
    
end
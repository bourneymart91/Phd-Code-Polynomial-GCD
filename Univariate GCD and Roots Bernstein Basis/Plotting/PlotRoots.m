global SETTINGS

% given the roots of f and g, plot them on a line
if(SETTINGS.PLOT_GRAPHS)
    
    figure_name = 'Exact roots of f(x), g(x) and d(x)';
    figure('name', figure_name)
    hold on
    title('Roots of f and g on the real interval')
    scatter(f_roots(:,1),ones(size(f_roots(:,1))),'s','DisplayName','Roots of f(x)')
    try
        scatter(g_roots(:,1),ones(size(g_roots(:,1))),'x','DisplayName','Roots of g(x)')
    catch
        fprintf('could not plot exact roots of g\n')
    end
    try
        scatter(d_roots(:,1),ones(size(d_roots(:,1))),'o','DisplayName','Roots of d(x)')
    catch
        fprintf('Could not plot exact roots of d.\n')
    end
    xlabel('Real')
    legend(gca,'show')
    hold off
    
end
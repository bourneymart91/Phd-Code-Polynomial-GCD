%% Plot Graphs

if(SETTINGS.PLOT_GRAPHS)
    
    figure('name','Residuals in SNTLN')
    hold on
    title('Residuals in SNTLN')
    xlabel('Iteration Number')
    ylabel('log_{10} Residuals')
    plot((1:1:ite),log10(condition),'-s','DisplayName','Residual')
    legend(gca,'show');
    hold off
    
    figure('name','Theta Variation over Newton Raphson Iterations')
    hold on
    title('Variation of \theta over Newton Raphson Iteration')
    xlabel('Iteration Number')
    plot((1:1:ite),log10(th1),'-s','DisplayName','\theta_{1}')
    plot((1:1:ite),log10(th2),'-s','DisplayName','\theta_{2}')
    legend(gca,'show');
    hold off
    
    figure('name','Alpha variation over Newton Raphson Iterations')
    hold on
    title('Alpha variation over Newton Raphson Iterations')
    xlabel('Iteration Number')
    plot((1:1:ite),log10(alpha),'-s','DisplayName','\alpha')
    legend(gca,'show');
    hold off
    
end


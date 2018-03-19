global SETTINGS
if( SETTINGS.PLOT_GRAPHS_RANK)
    
    % plot the minimum singular values
    figure_name = sprintf('%s : Minimum Singular Values of %s',mfilename,SETTINGS.SYLVESTER_MATRIX_TYPE);
    figure('name',figure_name)
    title('Minimum Singular Value for each subresultant matrix S_{k,k}')
    hold on
    plot(log10(vMinimumSingularValues),'-s','DisplayName','Preprocessed');
    %plot(log10(min_sing_val_vec_unproc),'-s','DisplayName','Unprocessed');
    xlabel('k : index of subresultant')
    legend(gca,'show')
    ylabel('log_{10} Minimum Singular Value')
    
    hold off
    
end
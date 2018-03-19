%% Graph Plotting
% Plot Graph of ratio of max min elements.


% Plot Graph of ratio of max : min element of the diagonal elements of R1 from the QR decompositions.
figure('name','GetDegree - MaxMin - Row Diags')
x = 1:min_mn;
plot(x,log10(vRatio_MaxMin_Diagonals_R),'red-s');
hold on
axis([1,min_mn,0,inf])
legend('Max:Min diag element of subresultant S_{k}');
title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
ylabel('log_{10} max:min diag element')
hold off


% Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
figure('name','GetDegree - MaxMin - Row Norms')
x = 1:min_mn;
plot(x,log10(vRatio_MaxMin_RowNorm_R),'red-s');
hold on
axis([1,min_mn,0,inf])
legend('Max:Min Row Norms of Rows in R1 from the QR decomposition of S_{k}');
title('Max:Min Row Norms of Rows in R1 from the QR Decomposition of S_{k} (Original)');
hold off


% Plot graph of norms of each row (N) from the qr decompostion of each S_{k}

figure('name','GetDegree - RowNorm')
plot(Data_RowNorm(:,1),(log10(Data_RowNorm(:,2))),'*')
axis([0.9,min_mn,-inf,+inf])
xlabel('k')
ylabel('Normalised Row Sums of R1 in S_{k}')
title(['Normalised Row Sums of R1 fom the QR decomposition of each subresultant S_{k} \newline '...
    'm = ' int2str(m) ', n = ' int2str(n) '(Original)']);
hold off


figure(4)
plot(Data_DiagNorm(:,1),(log10(Data_DiagNorm(:,2))),'*')
axis([0.9,min_mn,-inf,+inf])
xlabel('k')
ylabel('Normalised Diagonals of R1 in S_{k}')
title(['Normalised Diagonals in R1 matrix from the QR decomposition of each subresultant S_{k} \newline '...
    'm = ' int2str(m) ', n = ' int2str(n) '(Original)']);
hold off


% Plot graph of values of alpha for each subresultant

figure('name','GetDegree - Alphas')
plot(1:1:length(alpha_vec),log10(alpha_vec),'-s')
hold on
xlabel('k')
ylabel('log_{10} \alpha')
title('Optimal values of \alpha for each subresultant S_{k} (Original)')
hold off


% Plot graph of values of theta for each subresultant

figure('name','GetDegree - Thetas')
plot(1:1:length(theta_vec),log10(theta_vec),'-s')
hold on
xlabel('k')
ylabel('log_{10} \theta')
title('Optimal values of \theta for each subresultant S_{k} (Original)')
hold off

% Plot graph of Residuals by QR and SVD, if using the residual method to
% calculate the degree of the GCD.

x = 1:1:k;
figure(7)
plot(x,log10(vMinimumResidual),'red-o','DisplayName','Residuals by QR');
hold on
axis([1,min_mn,-inf,+inf])
ylabel('log_{10} Residual')
xlabel('k')
title('Residual obtained by removing optimal column from S_{k} (Original)');
legend(gca,'show')
hold off


switch PLOT_GRAPHS
    case 'y'
        % Plot for report
        figure('name','Get Degree - Max:Min - Row Norms of R1')
        hold on
        plot(log10(vRatio_MaxMin_RowNorm_R),'-s')
        xlabel('k : index of subresultant S')
        ylabel('log_{10} ratio max:min row sum r_{i} in R1')
        title('Plotting max:min row sum of the rows r_{i} of R1 from the QR decomposition of each subresultants S_{k}')
        hold off
        
        % Plot
        figure('name','Get Degree - Minimum Singular Values')
        hold on
        plot(log10(vMinimumSingularValues),'-s')
        hold off
        
        % Plot for report
        figure('name','Get Degree - Max:Min - Diagonals of R1')
        hold on
        plot(log10(vRatio_MaxMin_Diagonals_R),'-s')
        xlabel('k : index of subresultant S_{k}');
        ylabel('log_{10} ratio max:min diagon entries of R1')
        title('Plotting max:min diagonal entries of R1 from the QR decomposition of each subresultants S_{k}')
        hold off
    case 'n'
    otherwise
        error('error:')
end

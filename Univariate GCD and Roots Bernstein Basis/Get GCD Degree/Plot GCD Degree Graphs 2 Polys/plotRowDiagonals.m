function plotRowDiagonals(arr_RowDiagonals, limits_k, limits_t)
%
% % Inputs
%
% arr_RowDiagonals : (Array of Vectors) containing the diagonals of the matrices R_{k}, 
% from the QR decomposition of S_{k} for k = lower_lim:upper_lim 
%
% limits_k : [Int Int] : Limits on the possible values of k
%
% limits_t : [Int Int] : Limits 


%global SETTINGS

lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

% Plot graph of norms of each row (N) from the qr decompostion of each S_{k}
figure_name = sprintf('%s : Row Sum Norm',mfilename);
figure('name',figure_name)
hold on

for i = 1 : 1 : length(arr_RowDiagonals)
    
    % Get set of diagonals
    vRowDiags = arr_RowDiagonals{i};
    
    vec_i = i.*ones(length(vRowDiags),1);
    
    plot(vec_i, log10(vRowDiags) ,'*')

end

xlabel('k')
%ylabel(sprintf('Diagonals of R1 from QR decomposition of %s',SETTINGS.SYLVESTER_MATRIX_VARIANT))
%title(sprintf('Diagonals of R1 from QR decomposition of %s',SETTINGS.SYLVESTER_MATRIX_VARIANT));

grid on

xlim([lowerLimit_k upperLimit_k]);
%vline(limits_t,{'r','r'})


xlabel('$k$','Interpreter','latex','FontSize',20);
ylabel('$\log_{10} \left( R_{1,k}(i,i) \right)$','Interpreter','latex','FontSize',20);
hold off


% Figure size and location
myplot = gca;
myval_side = 0.10;
myval_base = 0.08;
set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])
set(gcf, 'Position', [100, 100, 710, 650])

box on
grid on
hold off
end
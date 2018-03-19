function plotMaxMinRowDiagonals(vMaxDiagR1, vMinDiagR1, limits_k, limits_t, rank_range)
% % Plot the max : min ratio of diagonals of matrices R_{1,k} for a range
% of k values, where k is the index of the k-th subresultant matrix and
% R_{1, k} is obtained by QR decomposition of S_{k}(f,g)
%
% % Inputs
%
% vMaxDiagR1 : (Vector) Maximum value of set of diagonals of the matrix 
% R_{1, k} obtained by the QR decomposition of S_{k}, for a set of k values
%
% vMinDiagR1 : (Vector)  Minimum value of set of diagonals of the matrix
% R_{1, k} obtained by the QR decomposition of S_{k}, for a set of k values
%
% limits_k : [Int Int] Range of k values for which the Condition numbers of
% S_{k} have been computed
%
% limits_t : [Int Int] Range of k which are candidates for the degree of
% the GCD of f(x) and g(x)
%
% rank_range : [Float Float]

% % Get lower and upper limit of k values.
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);


% 
x_vec = lowerLimit_k : 1 : upperLimit_k;

% Plot Graph of ratio of max : min element of the diagonal elements of R1 from the QR decompositions.
figure_name = sprintf('%s : Max:min Row Diagonals',mfilename);
figure('name',figure_name)

vRatio_MaxMin_Diagonals_R = vMinDiagR1./vMaxDiagR1;

plot(x_vec, log10(vRatio_MaxMin_Diagonals_R), 'red-s');
xlim([lowerLimit_k upperLimit_k]);

hline([rank_range mean(rank_range)],{'-r', '-r', '-b'})

vline(limits_t,{'-r', '-r'})

grid on
hold on
legend('Max:Min diag element of subresultant S_{k}');
title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
ylabel('log_{10} max:min diag element')
hold off

end




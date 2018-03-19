function [] = plotMaxMinDiagonalR1(vRatio_MaxMin_DiagonalEntry, limits_k, limits_t, rank_range)
%
% % Inputs
%
% vRatio_MaxMin_DiagonalEntry : (Vector)
%
% limits_k : (Int Int)
%
% limits_t : (Int Int)
%
% rank_range : [Float Float]

global SETTINGS

%
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);


figure_name = [mfilename ' : ' sprintf('Max;min Diagonals of R1 from QR Decomposition of %s', SETTINGS.SYLVESTER_MATRIX_VARIANT)];
    
figure('name', figure_name)
hold on
x_vec = lowerLimit_k : 1 : upperLimit_k;
plot(x_vec, log10(vRatio_MaxMin_DiagonalEntry),'-s')

%
hline([rank_range mean(rank_range)],{'-r', 'r', 'r'});
vline(limits_t, 'r', 'r');



hold off


end
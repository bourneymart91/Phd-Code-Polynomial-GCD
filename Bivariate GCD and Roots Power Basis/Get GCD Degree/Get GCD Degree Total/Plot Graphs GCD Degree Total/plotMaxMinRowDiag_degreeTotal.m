function [] = plotMaxMinRowDiag_degreeTotal(vRatio_MaxMin_Diags_R1, limits_k, limits_t, rank_range)
%
% % Inputs
%
% vRatio_MaxMin_Diags_R1
%
% limits_k : [Int Int] My limits on degree of GCD
%
% limits_t : [Int Int] Limits on degree of GCD
%
% rank_range : [Float Float]

myLowerLimit = limits_k(1);
myUpperLimit = limits_k(2);

figure_name = sprintf([mfilename ' : ' 'Plot Ratio Max over Min Diagonals of R1']);
figure('name',figure_name)
hold on
vec_x = myLowerLimit : 1 : myUpperLimit;
plot(vec_x, log10(vRatio_MaxMin_Diags_R1),'-s');


hline([rank_range mean(rank_range)],{'-r','-r','-r'});
vline(limits_t, {'-r','-r'});

hold off

end
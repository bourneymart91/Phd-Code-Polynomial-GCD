function plotMaxMinRowNorm_degreeTotal(vRatio_MaxMin, limits_k, limits_t, rank_range)
%
% % Inputs
%
% vRatio_MaxMin : [Vector]
% 
% limits_k : [Int Int]
%
% limits_t : [Int Int]
%
% rank_range : [Float Float]


lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

figure_name = sprintf([mfilename ' : ' 'Plot Max : Min Row Norm']);
figure('name',figure_name)
hold on
vec_x = lowerLimit_k : 1 : upperLimit_k;
plot(vec_x, log10(vRatio_MaxMin),'-s')

vline(limits_t, {'-r','-r'})
hline([rank_range mean(rank_range)],{'-r','-r','-r'});

hold off

end

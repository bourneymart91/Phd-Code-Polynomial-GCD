function [] = plotMaxMinRowNorms(vMaxRowNormR1, vMinRowNormR1, limits_k, limits_t, rank_range)
%
% % Inputs
%
% vMaxRowNormR1 : (Vector)
% 
% vMinRowNormR1 : (Vector)
%
% limits_k : (Int Int)
% 
% limits_t : (Int Int)
%
% rank_range : [Float Float]

% Get upper and lower limit of k
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);


vec_x = lowerLimit_k:1:upperLimit_k;

vMetric = vMinRowNormR1 ./ vMaxRowNormR1;


figure_name = '';
figure('name',figure_name)
hold on
plot(vec_x, log10(vMetric));

hline([rank_range mean(rank_range)], {'-r','-r','-b'})
vline(limits_t, {'-r','-r'});




hold off



end

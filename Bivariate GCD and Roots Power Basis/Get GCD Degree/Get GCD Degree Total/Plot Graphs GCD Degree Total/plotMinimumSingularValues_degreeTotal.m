function plotMinimumSingularValues_degreeTotal(vMinimumSingularValues, limits_k, limits_t, rank_range)
%
% % Inputs
%
% vMinimumSingularValues : (Vector) Vector of minimum singular values of
% each S_{k}
%
% limits_t : [(Int) (Int)]
%
% limits_t : [(Int) (Int)]
%
% rank_range : [Float Float]


lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

figure_name = sprintf([mfilename ' : ' 'Minimum Singular Values']);
figure('name',figure_name);
hold on

xlabel('k');
ylabel('log_{10} Singular Values');
x_vec = lowerLimit_k : 1 : upperLimit_k;
xlim(limits_k);
plot(x_vec, log10(vMinimumSingularValues), '-s');

hline([rank_range mean(rank_range)],{'-r','-r','-r'});
vline(limits_t, {'-r','-r'});

hold off

end
function plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t, rank_range)
%
% % Inputs 
% 
% vMinimumSingularValues : (Vector)
%
% limits_k : (Int) (Int) : Limits on the degree of the GCD. Defines the range
% of k values 
%
% limits_t : [Int Int] : Prior computed limits for the degree of the GCD
%
% rank_range : [Float Float] 
%


% 
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

%
lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

%
rank_range_low = rank_range(1);
rank_range_high = rank_range(2);

%
x_vec = lowerLimit_k : 1 : upperLimit_k;

figure_name = sprintf('%s : Minimum Singular Values of S_{k}',mfilename);
figure('name',figure_name);
hold on
plot(x_vec, log10(vMinimumSingularValues), '-s', 'DisplayName', 'Singular Values','LineWidth', 2)
ylabel('log \sigma(k)')
xlabel('k')

try
   xlim([1, length(x_vec)]) 
catch
end


% Plot vertical lines
vline(lowerLimit_t)
vline(upperLimit_t)

% plot horizontal lines
hline(rank_range_low);
hline(rank_range_high);


% Figure size and location
myplot = gca;
myval_side = 0.10;
myval_base = 0.08;
set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])
set(gcf, 'Position', [100, 100, 710, 650])

grid on
box on
xlabel('$k$','Interpreter', 'latex', 'FontSize', 20)
ylabel('$\log_{10} \left( \sigma_{k} \right) $','Interpreter', 'latex', 'FontSize', 20)
hold off

end
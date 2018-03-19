function plotMaxMinSingularValues(vMaximumSingularValues, vMinimumSingularValues, limits_k, limits_t, rank_range)
% Plot the minimum/maximum singular values
%
% % Inputs
%
% vMaximumSingularValues : (Vector) Vector of Maximum singular values from
% each S_{k}(f,g)
%
% vMinimumSingularValues : (Vector) Vector of Minimal singular values from
% each S_{k}(f,g)
%
% limits_k : [Int Int] The range of values over which the Sylvester
% subresultants are constructed
%
% limits_t : [Int Int] The range of values in which the degree 't' of the
% GCD must lie.
%
% rank_range : [Float Float]

%
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

global SETTINGS

figure_name = sprintf('%s : max:min Singular Values of %s',mfilename, ...
    SETTINGS.SYLVESTER_MATRIX_VARIANT);

figure('name',figure_name);
hold on
try
    xlim([lowerLimit_k upperLimit_k]);
catch
end
x_vec = lowerLimit_k : 1 : upperLimit_k;
vMetric = log10(vMinimumSingularValues) - log10(vMaximumSingularValues);
plot(x_vec,vMetric,'-s','DisplayName','Singular Values')



hline([rank_range mean(rank_range)],{'-r','-r','-b'})
vline(limits_t,{'r','r'})


% Plot Labels
ylabel('log (\sigma_{k})')
xlabel('k')
legend(gca,'show');

grid on

hold off

end

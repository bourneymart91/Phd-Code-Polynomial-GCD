function plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t, rank_range)
%
% % Inputs 
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

figure_name = sprintf('%s : Minimum Singular Values of %s',mfilename, ...
    SETTINGS.SYLVESTER_MATRIX_VARIANT);

figure('name',figure_name);
hold on
try
xlim([lowerLimit_k upperLimit_k]);
catch
end
x_vec = lowerLimit_k : 1 : upperLimit_k;
plot(x_vec,log10(vMinimumSingularValues),'-s','DisplayName','Singular Values')



hline([rank_range mean(rank_range)],{'-r','-r','-b'})
vline(limits_t,{'r','r'})


% Plot Labels
ylabel('$ \log_{10} \left( \sigma(k) \right)$','Interpreter','latex')
xlabel('$k$','Interpreter','latex')
legend(gca,'show');

grid on

hold off

end

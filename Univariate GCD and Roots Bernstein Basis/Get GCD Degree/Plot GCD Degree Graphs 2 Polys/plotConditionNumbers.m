function [] = plotConditionNumbers(vConditionNumbers, limits_k, limits_t)
% Plot the condition number of the set of subresultant matrices S_{k}(f,g)
% 
% % Inputs
%
% vConditionNumbers : (Vector of Floats) Condition number of the
% subresultant matrices S_{k}
%
% limits_k : [Int Int] Range of k values for which the Condition numbers of
% S_{k} have been computed
%
% limits_t : [Int Int] Range of k which are candidates for the degree of
% the GCD of f(x) and g(x)

% Get lower and upper limit of k values.
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

global SETTINGS

figure_name = sprintf('%s : Condition Numbers of %s', mfilename, ...
    SETTINGS.SYLVESTER_MATRIX_VARIANT);

figure('name',figure_name);
hold on
try
xlim([lowerLimit_k upperLimit_k]);
catch
end
x_vec = lowerLimit_k : 1 : upperLimit_k;


plot_name = strcat('$\kappa$ : BoolPreproc : ' , num2str(SETTINGS.BOOL_ALPHA_THETA));
plot(x_vec, log10(vConditionNumbers), '-s','DisplayName', plot_name)

vline(limits_t,{'r','r'})


% Plot Labels and legends
ylabel('$\log_{10} \left( \kappa_{k} \right)$', ...
    'Interpreter', 'latex', ...
    'FontSize',20)


xlabel('$k$', 'Interpreter', 'latex', 'FontSize',20)

l = legend(gca, 'show');
set(l,{'Interpreter', 'FontSize','Location'},{'latex', 20, 'southwest'});


% Appearance
grid on
box on


% Setting Size of plot within window
myplot = gca;
myval_side = 0.10;
myval_base = 0.08;
set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])


myWindow = gcf;
windowWidth = 700;
windowHeight = 610;
set(myWindow, 'Position', [ 100 100 windowWidth windowHeight])


hold off
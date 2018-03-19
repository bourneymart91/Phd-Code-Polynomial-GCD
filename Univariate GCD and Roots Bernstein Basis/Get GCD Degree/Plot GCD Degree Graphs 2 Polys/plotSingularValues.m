function plotSingularValues(arr_SingularValues, limits_k, limits_t)
%
% % Inputs
%
% arr_SingularValues
%
% limits_k : [Int Int] Upper and lower bound of k
%
% limits_t : [Int Int] Upper and lower bound of t
%
% rank_range : [Float Float]

% Get limits of k values
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

% Get number of Sylvester subresultant matrices
nSubresultants = upperLimit_k - lowerLimit_k + 1;

global SETTINGS

figure_name = sprintf('%s : All Singular Values of %s ',mfilename, SETTINGS.SYLVESTER_MATRIX_VARIANT);
figure('name', figure_name);
hold on

try
    xlim([lowerLimit_k upperLimit_k])
catch
    
end

for i = 1: 1: nSubresultants
   
    % Get k
    k = lowerLimit_k + (i-1);
    
    vSingularValues = arr_SingularValues{i};
    
    vec_k = k.*ones(length(vSingularValues),1);
    
    plot(vec_k, log10(vSingularValues), '*')
    
end



% Plot vertical lines
vline(limits_t,{'r','-r'});

grid on
box on

xlabel('$k$','Interpreter','latex','FontSize',20);
ylabel('$\log_{10} \left( \sigma_{k,i} \right)$','Interpreter','latex','FontSize',20);
hold off


% Figure size and location
myplot = gca;
myval_side = 0.10;
myval_base = 0.08;
set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])
set(gcf, 'Position', [100, 100, 710, 650])

box on
grid on
end
 function [] = plotSingularValues(arr_SingularValues, limits_k, limits_t)
%
% % Inputs
%
% arr_SingularValues : (Array of Vectors)
%
% limits_k : [Int Int]
%
% limits_t : [Int Int]

lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

nSubresultants = upperLimit_k - lowerLimit_k + 1;

global SETTINGS
figure_name = sprintf('Singular Values of %s' ,SETTINGS.SYLVESTER_MATRIX_VARIANT);
figure('name',figure_name)



hold on
for i = 1:1:nSubresultants
   
    k = lowerLimit_k + (i-1);
    
    vSingularValues = arr_SingularValues{i};
    
    nSingularValues = length(vSingularValues);
    
    vec_k = k.*ones(nSingularValues,1);
    
    plot(vec_k, log10(vSingularValues),'*');

end

vline(limits_t,{'r','r'});

% Font Size
set(gca, 'FontSize', 10)

% Labels
xlabel('$k$','Interpreter','latex','FontSize', 20)
ylabel('$\log_{10}\left(  \sigma_{k,i}  \right)$','Interpreter', 'latex', 'FontSize', 20)

% markers



% Resizing Figure
myplot = gca;
myval_side = 0.12;
myval_base = 0.10;
set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])

% Set window size
set(gcf, 'Position', [100, 100, 600, 600])

% Display Properties
grid on
box on




hold off




end
function [] = PlotCoefficients(arrPolys, arrNames)
%
% % Inputs
%
% arrPolys : (Array of Matrices)
%
% arrNames : (Array of Strings)


figure()
hold on
for i = 1 : 1 : length(arrPolys)
    
    
    fxy = arrPolys{i};
    m = GetDegree_Bivariate(fxy);
    
    nCoefficients_fxy = nchoosek(m + 2, 2);
    vFxy = GetAsVector(fxy);
    vFxy = vFxy(1 : nCoefficients_fxy);
    
    x_vec = 1 : 1 : nCoefficients_fxy;
    
    polyName = arrNames{i};
    plot(x_vec, log10(vFxy), 'DisplayName', polyName)
    
end

grid on
l = legend(gca,'show');
set(l,'Interpreter', 'latex');
set(l,'FontSize',20);
hold off



% Labels
xlabel('$i$','Interpreter','latex','FontSize', 20)
ylabel('$\log_{10}\left(  \Re  \right)$','Interpreter', 'latex', 'FontSize', 20)


% Resizing Figure
myplot = gca;
myval_side = 0.12;
myval_base = 0.10;
set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])

% Set window size
set(gcf, 'Position', [100, 100, 600, 600])


end
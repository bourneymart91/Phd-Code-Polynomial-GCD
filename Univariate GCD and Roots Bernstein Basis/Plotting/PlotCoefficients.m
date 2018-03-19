
function [] = PlotCoefficients(arrPolys, arrLabels, arrStyles)
% Plot coefficients of the set of polynomials stored in the array 'arrPolys'
%
% % Inputs
%
% arrPolys : (Array of vectors) Array of coefficients of polynomials
%
% arrLabels : (Array of Strings) Array of labels for each polynomial
%
% arrStyles : (Array of Strings) Style of each of the polynomials


global SETTINGS
if SETTINGS.PLOT_GRAPHS
    
    % Plot on log scale
    bool_log = true;
    
    
    % Get number of polynomials in the array
    nPolys = length(arrPolys);
    
    % Initalise a figure
    figure()
    hold on
    
    % Initliase a vector to store the degree of each polynomial
    vDegree = zeros(nPolys,1);
    
    
    % For each polynomial, plot its coefficients
    for i = 1 : 1 : nPolys
        
        myStyle = arrStyles{i};
        
        % Get vector of coefficients of the i-th polynomial
        fx = arrPolys{i};
        
        % Get name of the i-th polynomial
        name = arrLabels{i};
        
        % Get the degree of the i-th polynomial
        vDegree(i) = GetDegree(fx);
        
        % Initialise vector of x values
        vec_x = 0 : 1 : vDegree(i);
        
        
        if (bool_log == true)
            fx = log10(abs(fx));
        end
        
        % Plot i-th polynomial
        plot(vec_x, fx, myStyle, 'DisplayName',name,'LineWidth',2)
        
    end
    
    xlim([1, max(vDegree)]);
    
    % Labels and legends set up
    xlabel('$i$ : Coefficient Index','Interpreter','latex', 'FontSize', 20);
    ylabel('$\log_{10} \left( \Re \right)$', 'Interpreter', 'latex', 'FontSize',20);
    l = legend(gca,'show');
    set(l,{'Interpreter','FontSize','Location'},{'latex',20, 'southwest'});
    hold off
    
    
    % Figure size and location set up
    myplot = gca;
    myval_side = 0.10;
    myval_base = 0.08;
    set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])
    set(gcf, 'Position', [100, 100, 710, 650])
    
    box on
    grid on
    
end



end

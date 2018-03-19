function [] = Plot_fx(fx, a, b, figure_name)
% Plot the Bernstein polynomial f(x) and its control points defined over
% the interval [a,b]
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x) 
%
% a : (Float) Lower end of interval
%
% b : (Float) Upper end of interval
%
% figure_name : (String)


global SETTINGS

switch SETTINGS.PLOT_GRAPHS
    case true
        
        % Get the control points of f(x)
        Pk = GetControlPoints(a, b, fx);
        
        % Initialise the column vector of x values
        x = linspace(a, b, 100)';
        
        % Get number of entries in x vector
        nEntries_x = size(x, 1);
        
        % Initialise a y vector of the same size as x
        y = zeros(nEntries_x, 1);
        
        % Set the y ordinate values.
        for i = 1:1:nEntries_x
            
            % Evaluate f(x) at x_{i}
            y(i) = Bernstein_Evaluate(fx,x(i));
            
        end
        
        
        figure('name',figure_name)
        plot(x,y,'-');
        hold on
        grid on
        scatter(Pk(:,1),Pk(:,2));
        hold off
        
    case false
        
    otherwise
        error('err');
end
end
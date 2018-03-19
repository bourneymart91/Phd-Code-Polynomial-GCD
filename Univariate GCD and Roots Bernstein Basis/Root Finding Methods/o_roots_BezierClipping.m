function [root_mult_array] = o_roots_BezierClipping(fx)
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% % Outputs
%
% root_mult_array : [Float Int] Matrix consisting of two columns containing
%   (i) a root and (ii) its multiplicity.



% Pseudocode from 'Computing roots of polynomials by quadratic clipping'
% Barton and Juttler
%
% Algorithm bezclip(p,[a,b])
% where p polynomial
% 1. if length of interval [a,b] > e then
% 2.    C <- Convex hull of control points of p with respect to [a,b]
% 3.    if C intersects t axis then
% 4.        Find [a',b'] by intersecting C with the t axis
% 5.        if |a'-b'|< 1/2*|a-b| then
% 6.            return (bezclip(p,[a',b'])
% 7.        else
% 8.            return (bezclip(p,[a,1/2(a+b)]) union bezclip(p,1/2(a+b,b))
% 9.        end if
% 10.   else
% 11.   return empty set
% 12.   end if
% 13. else
% 14. return ([a,b])
% 15. end if


% bool_preproc = '';
% low_rank_approx_method = '';
% apf_method = '';
% SetGlobalVariables(bool_preproc,low_rank_approx_method,apf_method);

% Set interval size
global xmin
global xmax
global ymin
global ymax


interval_low = 0;
interval_high = 1;

% max error
epsilon = 1e-5;

% Set iteration number to one
outer_loop_ite = 1;

% Initialise a vector to store calculated roots
root_mult_array = [];

% While f(x) is not a constant
while GetDegree(fx) >= 1
    
    if GetDegree(fx) == 1
        

        % Curve is a straight line.
        m = 1;
        x = linspace(interval_low, interval_high, m + 1)';
        CP = [x, fx];
        
        xmin = min(CP(:,1));
        xmax = max(CP(:,1));
        ymin = min(CP(:,2));
        ymax = max(CP(:,2));
        
        x_intercept_new = GetXIntercept(CP);
        fprintf('x intercept : %2.4f', x_intercept_new)
    else
        
        % Get degree of polynomial f(x)
        m = GetDegree(fx);
        
        % Define the set of control points of f(x)
        % Split the unit interval into m+1 parts
        x = linspace(interval_low, interval_high, m + 1)';
        CP = [x,fx];
        
        xmin = min(CP(:, 1));
        xmax = max(CP(:, 1));
        ymin = min(CP(:, 2));
        ymax = max(CP(:, 2));
        
        
        x_intercept_old = interval_low;
        x_intercept_new = interval_high;
        
        
        inner_loop_ite = 1;
        
        while abs(x_intercept_new - x_intercept_old) > epsilon
            
            fprintf('New Intercept : %f \n', x_intercept_new)
            fprintf('Old Intercept : %f \n', x_intercept_old)
            
            fprintf('Iteration %i \n',inner_loop_ite)
            
            % Get the convex hull of the control points of f(x)
            k = convhull(CP(:,1),CP(:,2));
            
            
            fprintf('Convex hull over interval %f - %f \n',CP(1,1), CP(end,1));
            
            convex_hull_vertices = flipud(CP(k,:));
            
            % for each line of the convex hull, check if it intersects
            % Plot the control points
            %figure_name = sprintf('%s : Plotting f(x)',mfilename);
            %Plot_fx(CP,0,1,'')
            
            % Get the first point at which the convex hull crosses the x axis.
            x_intercept_old = x_intercept_new;
            
            x_intercept_new = GetXIntercept(convex_hull_vertices);
           
            if(x_intercept_new == -1000)
                break;
            end
            
            % Evaluate f(x) at the point.
            Bernstein_Evaluate(fx, x_intercept_new);
            
            % Subdivide at the point t1.
            [~, Pk_right] = BezierSubdivide(CP, m, x_intercept_new);
            
            CP = Pk_right;
            inner_loop_ite = inner_loop_ite + 1;
        end
    end
    % Save all plots
    
    if (x_intercept_new == -1000)
        break;
    end
    
    %dir_name = sprintf('outputs-%i',outer_loop_ite);
    %save_all_figures_to_directory(dir_name);
    
    % Close all open plots
    %close all;
    
    fprintf('** Root found at : %f',x_intercept_new)
    
    % Get the root.
    root = x_intercept_new;
    
    % Add the root to the root array.
    root_mult_array = [root_mult_array ; root 1 ];
    
    % Get polynomial of the factor (x-r)
    tx =[...
        -root;
        1-root;
        ];
    
    % Deconvolve root from polynomial
    fx = Deconvolve(fx, tx);
    Plot_fx(fx,0,1,'')
    
    % Increment iteration number
    outer_loop_ite = outer_loop_ite + 1;
    
end


end






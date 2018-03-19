function [] = o_Intersection_2DCurveExplicit_2DCurveExplicit(...
    ex_num_f, ex_num_g, el, bool_preproc, low_rank_approx_method)
% o_Intersection_2DCurveExplicit_2DCurveExplicit...
%       (ex_num_f,ex_num_g,bool_preproc,low_rank_approx_method)
%
% Given two explicitly defined curves f(x) and g(x) in power form,
% calculate their intersections.
%

SetGlobalVariables(bool_preproc, low_rank_approx_method);
global SETTINGS



% Get the type of example used
EXAMPLE_TYPE = 'Coefficients';

switch EXAMPLE_TYPE
    case 'Coefficients'
        switch ex_num_f
           
            case '1'
                % f = (x-2)^2 (x-3)^3 (2x+3)
                % f = 2 x^6-23 x^5+95 x^4-141 x^3-81 x^2+432 x-324
                fx = [-324; 432; -81; -141; 95; -23; 2];
                
                % g = (x-2)^2 (x-3)^3 (x+4)
                % g = x^6-9 x^5+15 x^4+97 x^3-468 x^2+756 x-432
                gx = [-432; 756; -468; 97; 15; -9; 1];
                
                % f-g = (x-2)^2 (x-3)^3 (x-1)
                % f-g = x^6-14 x^5+80 x^4-238 x^3+387 x^2-324 x+108
                % [108; -324; 387; -238;80;-14;1]


        end
    case 'Roots'
        
        % Get the roots of two implicitly defined polynomial curves.
        fx_root_mult_arr = Examples_Univariate_Implicit(ex_num_f);
        gx_root_mult_arr = Examples_Univariate_Implicit(ex_num_g);
        
        % Print the factorisation of f(x) and g(x)
        PrintFactorization(fx_root_mult_arr, 'f')
        PrintFactorization(gx_root_mult_arr, 'g')
               
        % Get the coefficients of the polynomails f(x) and g(x)
        fx = GetCoefficientsFromRoots(fx_root_mult_arr);
        gx = GetCoefficientsFromRoots(gx_root_mult_arr);
        
        % Print the coefficients of f(x) and g(x)             
        PrintCoefficientsBivariate(fx, 'f');
        PrintCoefficientsBivariate(gx, 'g');
        
        % Add noise to the coefficients
        fx = AddNoiseToPoly(fx, el);
        gx = AddNoiseToPoly(gx, el);
                
    case 'SamePoly'
        fx = Examples_Univariate_Implicit(ex_num_f);
        % Add noise to fx
        el = 1e-14;
        gx = fx + (el.*ones(size(fx)));
        PrintCoefficientsBivariate(fx,'f')
        PrintCoefficientsBivariate(gx,'g')
        
end

% Get degree of polynomial f(x)
m = GetDegree(fx);

% Get degree of polynomial g(x)
n = GetDegree(gx);

%%
%% Plot the explicitly defined curves.
t = -5:0.001:5;

f_y = polyval(flipud(fx),t);
g_y = polyval(flipud(gx),t);


IF(SETTINGS.PLOT_GRAPHS)
    
        figure('name','Curve Plot : f(y) and g(y)')
        hold on
        plot(t,f_y,'DisplayName','f(y)')
        plot(t,g_y,'DisplayName','g(y)')
        legend(gca,'show');
        hold off

    
end


%%
INTERSECTION_METHOD = 'Common Factors First';
INTERSECTION_METHOD = 'All Roots';

switch INTERSECTION_METHOD
    case 'All Roots'
        % Get the polynomial h(x) = f(x) - g(x)
        % For this, f(x) and g(x) must be of the same length.
        
        % Create a zeros vector of the length of h(x)
        zero_vector = zeros(min(m,n),1);
        
        % Replace f(x) by f(x) padded with zeros for coefficients
        % x^{m+1},...
        f = zero_vector;
        f(1:m+1) = fx;
        fx_inflated = f;
        
        % Replace g(x) by g(x) padded with zeros for coefficients
        % x^{n+1},...
        g = zero_vector;
        g(1:n+1) = gx;
        gx_inflated = g;
        
        % Get the coefficients h(x)
        hx = fx_inflated - gx_inflated;
        
        % The degree of fx - gx may be reduced, find the last non-zero
        % coefficient
        i2 = find(hx, 1, 'last'); 
        hx = hx(1:i2);

        
        % Check for zero coefficients
        if hx(end) == 0
            hx(end) = []
        end
        
        %fprintf('Roots by Matlab Method \n')
        %roots(flipud(hx))
        
    case 'Common Factors First'
        
        % Get the GCD of f(x) and g(x)
        [~,~,dx, ux, vx , alpha,theta, ~, lambda, mu] = o_gcd_mymethod(fx,gx);
        
        % Get the roots of polynomial d(x)
        [t,~] = size(dx);
        if t >1 
            my_roots = o_roots_mymethod(dx);
        end
        
        % Get degree of u(x)
        m = GetDegree(ux);
        
        % Get degree of v(x)
        n = GetDegree(vx);
        
        % % Subtract u(x) from v(x)
        
        vZeros = zeros(max(m,n),1);
        u = vZeros;
        u(1:m+1) = ux;
        
        v = vZeros;
        v(1:n+1) = vx;
        
        % Get the coefficients h(x).
        hx = u - v;
        
        % Check for zero coefficients
        if hx(end) == 0
            hx(end) = [];
        end
end

% Get the degree of h(x)
m = GetDegree(hx);
if m ~=0
    % Get the roots of h(x)
    roots2 = o_roots_mymethod(hx);
end


o_roots_matlab(fx);


%%

try
display(roots2)
catch
    fprintf('err')
end

% % Given the roots in x, calculat the corresonding y values

% Get number of roots
[nEntries,~] = size(roots2);

% Initialise a matrix to store (x,y) pairs
xy_pairs = zeros(nEntries,2);

% For each root
for i = 1:1:nEntries
    
    % Get x coordinate
    x_val = roots2(i);
    
    % Evaluate polynomial f(x) at x_{i}
    y_val = polyval(flipud(fx),x_val);
    
    % Add to matrix of intersection coordinate pairs.
    xy_pairs(i,:) = [x_val, y_val];
end

display(xy_pairs)



% Get the points of intersection.
% for each root in roots




end

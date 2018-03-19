function [] = o_Intersection_PlanarParametricCurve_PlanarParametricCurve(bool_preproc,low_rank_approx_method)
% Given two integral parametrically defined curves, obtain the points of
% intersection.

SetGlobalVariables()

% fxt : vector of coefficients of x(t) with increasing powers of t.
% fyt : vector of coefficients of y(t) with increasing powers of t.
C1_xt = [0.1; 0.8; -0.1];
C1_yt = [0.5; 3; -2.5];

C2_xt = [0; 2; 0.1];
C2_yt = [0; 1; 0.1];

% Get the implicit representation of the Integral parametric polynomial
% C_{1}(x(t),y(t))
C1_implicit = Implicitize_Integral_Parametric(C1_xt,C1_yt);

% Get implicit Representation, C1 is given by f(x,y) = 0.
fprintf('The Implicit Equation of C_{1}(x,y) is given by')
fprintf('\n')
disp(C1_implicit)


PrintCoefficientsBivariate(C1_implicit,'C1')

%%
% Get the number of rows and columns in C1(x,y)
[rowsC1,colsC1] = size(C1_implicit);

% Produce the polynomial C3 given by the susbstitution of x(t) and y(t) 
% into f(x,y) = 0. Where C3 is a curve with variable t.
C3 = 0;

% for each row in C_{1}(x,y)
for i = 0:1:rowsC1-1
    % for each column in C_{1}(x,y)
    for j = 0:1:colsC1-1
        
        x_component = 1;
        for k = 0:1:i
            x_component = conv(C2_xt,x_component);
        end
        
        y_component = 1;
        for k = 0:1:j
            y_component = conv(C2_yt,y_component);
        end
        
        uij = conv(x_component,y_component);
        
        uij = uij .* C1_implicit(i+1,j+1);
        
        C3 = PolyAdd(C3, uij);
    end
end

% Print the coefficients of the curve C3
fprintf('The Curve C_{3} is given by \n')
fprintf('\n')
PrintCoefficientsBivariate(C3,'C3')

% Get roots of C_{3} in terms of t
fprintf('roots by my method \n')
[root_mult_array] = o_roots_mymethod(C3);
display(root_mult_array);

% Get roots of C_{3} by my method
vRoots = roots(fliplr(C3'));

fprintf('The roots in terms of t are given by')
fprintf('\n')
disp(vRoots)


[nRoots,~] = size(vRoots);

x = zeros(1,nRoots);
y = zeros(1,nRoots);

% For each root r_{i} in the vector vRoots
for i = 0:1:nRoots-1
    
    % Get the ith root
    rt = vRoots(i+1);
    
    % Substitute in to x(t)
    x_sum = 0;
    
    % for each coefficient in C1_x(t)
    for j = 0:1:size(C2_xt,1)-1
        x_sum = x_sum + (C2_xt(j+1) .* (rt.^j));
    end
    
    x(i+1) = x_sum;
    
    y_sum = 0;
    for j = 0:1:size(C2_yt) -1
        y_sum = y_sum + (C2_yt(j+1) .* (rt.^j));
    end
    
    y(i+1) = y_sum;
    
    
end


fprintf('The intersections are given by [x,y]')
fprintf('\n')
disp([x' y'])


% Substitute root values back into C_{1} x(t) and y(t) to obtain set of
% intersection points
global SETTINGS
if(SETTINGS.PLOT_GRAPHS)

        

        figure('name','Plotting C_{1} and C_{2}')
        hold on
        
        % Get symbolic expressions for C1 = (x(t),y(t))
        x1t = poly2sym(flipud(C1_xt));
        y1t = poly2sym(flipud(C1_yt));
        
        % Get symbolic expressions for C2 = (x(t),y(t))
        x2t = poly2sym(flipud(C2_xt));
        y2t = poly2sym(flipud(C2_yt));
        
        % Plot the Symbolic Curves C1 and C2
        ezplot(x1t,y1t);
        ezplot(x2t,y2t);
        
        axis([-10,10,-10,10]);
        grid on
        hold off
        
   
end


end
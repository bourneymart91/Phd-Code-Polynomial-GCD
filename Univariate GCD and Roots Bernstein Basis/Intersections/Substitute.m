function [coef_Bernstein_Poly]= Substitute(CP_f, gxy)
% Substitute x(t) and y(t) from f(x,y) into g(x,y)

% Get degree of input polynomial 
[~,c] = size(CP_f);
n = c-1;

% Get coefficients of f(x)
f_x =  CP_f(1,:)';

% Get Coefficients of f(y)
f_y =  CP_f(2,:);

% Get Degree of g(x,y)
[n1,n2] = GetDegree_Bivariate(gxy);

% Initialise the univariate Bernstein basis polynomial obtained from the
% substitution.
coef_Bernstein_Poly = zeros((n^2)+1,1);

% for each row in the implicit representation of g(x,y)
for i = 0:1:n1
    
    % for each column of g(x,y)
    for j = 0:1:n2
        
        
        % Get the coefficient b_{i,j} in g(x,y)
        coef = gxy(i+1,j+1);
        
        % Coefficient b_{i,j} has x^{i}
        
        % Get x(t)^{i}
        x_component = 1;
        for k = 1:1:i   
            x_component = Bernstein_Multiply(x_component,f_x);
        end
        
        % Get y(t)^{j}
        y_component = 1;
        for k = 1:1:j
            y_component = Bernstein_Multiply(y_component',f_y')';
        end
        
        % Get x(t)^i * y(t)^j
        xy_comp = Bernstein_Multiply(x_component,y_component);
        
        
        uij =  xy_comp;
        
        % output polynomial will be of degree 2n
        m = n^2;
        
        % degree elevate uij
        [r1,~] = size(uij);
        curr_deg_uij = r1-1;
        num_deg_elv_req =  m-curr_deg_uij;

        uij = Bernstein_DegreeElevate_Univariate(uij,num_deg_elv_req);
        
        uij = coef .* uij;
        
        coef_Bernstein_Poly = coef_Bernstein_Poly + uij;
    end
end
end
function [h] = AddPolynomials(fx, gx)
% Add two polynomials f(x) and g(x)
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% gx : (Vector) Coefficients of polynomial g(x)
%
% % Outputs
%
% hx : (Vector) Coefficients of polynomial h(x)


% if nCols > 1, not a column vector
if (size(fx,2) >1 || size(gx,2) >1)
   error('Not a column vector') 
end

% Get degree of polynomial f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);


if m < n
    fx = [fx; zeros(n-m,1)];
elseif n < m
    gx = [gx; zeros(m-n,1)];
end

h = fx + gx;

end
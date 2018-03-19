function [m1,m2] = GetDegree_Bivariate(fxy)
% GetDegree_Bivariate(fxy)
%
% Get the degree of the bivariate polynomaial f(x,y).
%
% Inputs.
% 
% fxy : Matrix of coefficients of polynomial f(x,y)
%
% Outputs.
%
% m1 : Degree of polynomial f(x,y) with respect to x.
%
% m2 : Degree of polynomial f(x,y) with respect to y.

% Get dimensions of the matrix of coefficients of f(x,y).
[nRows,nCols] = size(fxy);

% Set m1 the degree with respect to x
m1 = nRows - 1;

% Set m2 the degree with respect to y
m2 = nCols - 1;

end
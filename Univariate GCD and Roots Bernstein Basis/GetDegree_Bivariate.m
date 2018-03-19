function [m1, m2] = GetDegree_Bivariate(fxy)
% GetDegree(fxy)
%
% Get the degree of the polynomial f(x,y), whose coefficients are given as 
% a matrix.
%
% % Inputs
%
%
% fxy : (Matrix) Matrix of coefficients of the polynomial f(x,y)
%
%
% % Outputs
%
%
% m1 : (Int) Degree of f(x,y) with respect to x
%
% m2 : (Int) Degree of f(x,y) with respect to y





% Get number of rows and number of columns in matrix of coefficients of f(x,y)
[nRows,nCols] = size(fxy);

% Get degree m1 of f(x,y) with respect to x
m1 = nRows - 1;

% Get degree m2 of f(x,y) with respect to y
m2 = nCols - 1;

end
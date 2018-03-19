function [m1, m2] = GetDegree_Bivariate(fxy)
% Get the degree of the bivariate polynomial f(x,y). Note that for all
% polynomials in the Bivariate triangular Bernstein form the degree m = m1
% = m2.
%
% Inputs
% 
% fxy : (Matrix) Coefficients of polynomial f(x,y).
%
% Outputs
%
% m1 : (Int) Degree of f(x,y) with respect to variable x.
%
% m2 : (Int) Degree of f(x,y) with respect to variable y.




% Get the number of rows and columns of the matrix containing the
% coefficients of f(x,y)
[nRows, nCols] = size(fxy);

m1 = nRows - 1;
m2 = nCols - 1;



end
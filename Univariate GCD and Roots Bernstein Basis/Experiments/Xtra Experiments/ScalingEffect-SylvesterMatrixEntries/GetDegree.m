function [m1, m2] = GetDegree(fxy)
% Get the degree of f(x,y) with respect to both x and y
%
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)

[r,c] = size(fxy);

% Get degree with respect to x
m1 = r - 1;

% Get degree with respect to y
m2 = c - 1;


end
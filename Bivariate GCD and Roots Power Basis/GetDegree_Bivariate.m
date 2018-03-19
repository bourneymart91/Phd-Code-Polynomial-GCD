function [m1, m2] = GetDegree_Bivariate(fxy)
% Get the degree of the polynomial f(x,y)

[r,c] = size(fxy);
m1 = r - 1;
m2 = c - 1;


end
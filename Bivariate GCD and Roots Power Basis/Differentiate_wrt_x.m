function [partial_fxy] = Differentiate_wrt_x(fxy, m)
% Given the polynomial fxy in matrix form
%
%           1   y   y^{2}
%          _           _
%       1 |             |
%       x |             |
%   x^{2} |_           _|
%
% Obtain the partial derivative with respect to x
%

% Get the degree of fxy with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);


% Create a vector to multiply coefficients by
mult_vec = (0:1:m1);

% Multiply the rows of fxy by the multiplication vector
partial_fxy = diag(mult_vec) * fxy ;

partial_fxy(1,:) = [];

if (m2 == m)
   partial_fxy(:,end) = [] ;
end

end
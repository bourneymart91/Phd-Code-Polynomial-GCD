function [symbolic_exp] = GetSymbolicFromCoefficients(fxy, var_x, var_y)
% Given the matrix of coefficients, get the symbolic representation of the
% polynomial f(x,y).
%
% Inputs.
%
% fxy_matrix : Matrix of coefficients of the polynomial f(x,y)
%
% var_x : String of name of variable eg 'x'
%
% var_y : String of name of variable eg 'y'

% Note that the var_x and var_y can be changed to change variables, for
% example this fucntion can be used to write the symbolic expression of a
% polynomial f(s,t) by GetSymbolicFromCoefficients(f,'s','t')

% Initialise symbols x and y
x = sym(var_x);
y = sym(var_y);


% Get the degree of the polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% % 
% 

mat_x = diag(x.^(0:1:m1));
mat_y = diag(y.^(0:1:m2));

symbolic_mat = mat_x * fxy * mat_y;

symbolic_exp = sum(sum(symbolic_mat));

end



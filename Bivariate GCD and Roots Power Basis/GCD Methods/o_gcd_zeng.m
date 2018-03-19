function [dxy_mat,uxy_mat,vxy_mat] = o_gcd_zeng(fxy,gxy,tol)
% Prepare the polynomials f(x,y) and g(x,y) for input into Zengs method.
% Note that the inputs must be polynomials in string form

addpath('NAClabForMatlab2013')

x = sym('x');
y = sym('y');

f = GetSymbolicFromCoefficients(fxy,x,y);
g = GetSymbolicFromCoefficients(gxy,x,y);

f = char(expand(f));
g = char(expand(g));


%f = '-45*x*y - 15*x^3*y - 20*x*y^2 + 27*x*y^3 + 9*x^3*y^3 + 12*x*y^4';
%g = '45*x^2*y^2 + 15*x^2*y^3 - 27*x^2*y^4 - 9*x^2*y^5';

[dxy_str,uxy_str,vxy_str] = PolynomialGCD(f,g,tol);

dxy_mat_zeng = PolynString2CoefMat(dxy_str,{'x','y'});
dxy_mat = GetCoefficientMatrix(dxy_mat_zeng);

uxy_mat_zeng = PolynString2CoefMat(uxy_str,{'x','y'});
uxy_mat = GetCoefficientMatrix(uxy_mat_zeng);

vxy_mat_zeng = PolynString2CoefMat(vxy_str,{'x','y'});
vxy_mat = GetCoefficientMatrix(vxy_mat_zeng);

end

function my_mat = GetCoefficientMatrix(coef_mat)
% Given the coefficient matrix is zeng format, convert to coefficient
% matrix in my format

% Zeng format
% Each column contains three values [x power; y power; coeff]

% My format:
% Matrix consists of coefficients, position of coefficient indicates
% power, eg position (2,3) has power x^2 y^3.

% Get number of columns of coef_mat
[~, nCols] = size(coef_mat);

% Get highest degree with respect to x
max_x_pwr = max(coef_mat(1,:));
max_y_pwr = max(coef_mat(2,:));

% Initialise my coefficient matrix
my_mat = zeros(max_x_pwr+1 , max_y_pwr +1);

% For each coefficient in Zengs matrix, put into my matrix
for i = 1:1:nCols
    x_pos = coef_mat(1,i) +1;
    y_pos = coef_mat(2,i) +1;
    my_mat(x_pos,y_pos) = coef_mat(3,i);
end
end
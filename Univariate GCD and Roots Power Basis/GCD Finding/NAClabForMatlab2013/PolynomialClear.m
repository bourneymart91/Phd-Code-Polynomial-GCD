%
% Clear tiny coefficients, tiny real parts, and tiny imaginary parts
% from a polynomial
%
% Syntax:  (pclear is the shortened alias of PolynomialClear)
%          >> g = PolynomialClear(f)
%          >> g = PolynomialClear(f,tol)
%
%   Input:   f --- (string or numeric) polynomial
%          tol --- (numeric)           threshold for being tiny
%
%  Output:   g --- (string or numeric) cleared polynomial
%
%  Example:
%
%   >> PolynomialClear('(3+2e-13*i) + 2.1e-15*x*y+(1e-14+2*i)*x^3',1e-10)
%
%   ans =
%
%   3 + (0+2i)*x^3

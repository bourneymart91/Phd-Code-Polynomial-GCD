function [px] = Bernstein_Multiply(fx, gx)
% Given the Coefficients of two polynomials f(x) and g(x) in vector form, 
% output the coefficients of the product p(x).
%
% % Input.
%
% fx : (Vector) Column vector of coefficients of the Bernstein Polynomial f(x)
%
% gx : (Vector) Column vector of coefficients of the Bernstein Polynomial g(x)
%
% % Output.
%
% px : (Vector) Column vector of coefficients of the Bernstein Polynomial p(x)



% Get the degree of polynomial g(x)
n = GetDegree(gx);

% Build the matrix D^{-1}*T(f)*Q1
DTQ = BuildDT1Q1(fx, n);

% Get the vector of coefficients of p(x) from the product f(x) and g(x).
px = DTQ*gx;

end
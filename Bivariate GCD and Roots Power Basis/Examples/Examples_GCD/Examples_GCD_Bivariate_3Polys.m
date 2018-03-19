function [fxy, gxy, hxy, uxy, vxy, wxy, dxy,...
    m, m1, m2,...
    n, n1, n2,...
    o, o1, o2,...
    t, t1, t2]  = Examples_GCD_Bivariate_3Polys(ex_num)
% Examples_GCD(ex_num)
% Get the coefficient matrix of the polynomials f(x,y), g(x,y), the GCD
% d(x,y), and cofactors u(x,y) and v(x,y) as well as their corresponding
% degrees.
%
% This file either produces an example from coefficient matrices or 
% 'from roots' where the coefficient matrices are generated.
%
% % Inputs
% 
% ex_num : Example number
%
% % Outputs
%
% fxy : (Matrix) Coefficient matrix of polynomials f(x,y)
%
% gxy : (Matrix) Coefficient matrix of polynomials g(x,y)
%
% hxy : (Matrix) Coefficient matrix of polynomials h(x,y)
%
% uxy : (Matrix) Coefficient matrix of polynomials u(x,y)
% 
% vxy : (Matrix) Coefficient matrix of polynomials v(x,y)
%
% wxy : (Matrix) Coefficient matrix of cofactor polynomials w(x,y)
%
% dxy : (Matrix) Coefficient matrix of polynomial d(x,y)
%
% m, m1, m2 : (Int) (Int) (Int) : Degree structure of polynomial f(x,y)
%
% n, n1, n2 : (Int) (Int) (Int) : Degree structure of polynomial g(x,y)
%
% o, o1, o2 : (Int) (Int) (Int) : Degree structure of polynomial h(x,y)
%
% t, t1, t2 : (Int) (Int) (Int) : Degree structure of polynomial d(x,y)

EXAMPLE_TYPE = 'From Coefficients';

switch EXAMPLE_TYPE

    case 'From Coefficients'
        [fxy, gxy, hxy, uxy, vxy, wxy, dxy,...
            m, m1, m2,...
            n, n1, n2,...
            o, o1, o2,...
            t, t1, t2] = Examples_GCD_FromCoefficients_Bivariate_3Polys(ex_num);
        
end
end



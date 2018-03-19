function [fxy, gxy, uxy, vxy, dxy,...
    m, m1, m2, n, n1, n2, t, t1, t2]  = Examples_GCD_Bivariate_2Polys(ex_num)
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
% ex_num : (String) Example Number
%
% % Outputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% uxy : (Matrix) Coefficients of cofactor polynomial u(x,y) 
%
% vxy : (Matrix) Coefficients of cofactor polynomial v(x,y) 
%
% dxy : (Matrix) Coefficients of GCD d(x,y)
%
% m m1 m2 : (Int) (Int) (Int) : Degree structure of f(x,y)
%
% n n1 n2 : (Int) (Int) (Int) : Degree structure of g(x,y)
%
% t t1 t2 : (Int) (Int) (Int) : Degree structure of d(x,y)

EXAMPLE_TYPE = 'From Coefficients';

switch EXAMPLE_TYPE
    case 'From Roots'
        
        % Get example polynomials
        [fxy, gxy,...
            uxy,vxy,...
            dxy,...
            m,m1,m2,...
            n,n1,n2,...
            t,t1,t2] = Examples_GCD_FromRoots_Bivariate(ex_num);
        
    case 'From Coefficients'
        [fxy, gxy,...
            uxy,vxy,...
            dxy,...
            m,m1,m2,...
            n,n1,n2,...
            t,t1,t2] = Examples_GCD_FromCoefficients_Bivariate_2Polys(ex_num);
        
end
end



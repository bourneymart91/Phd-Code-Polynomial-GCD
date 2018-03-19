function [fxy, m] = Examples_Roots_Bivariate(ex_num)
% Given an example number, get a set of roots, Build the polynomial and
% output the matrix of coefficients.
%
% % Inputs.
%
% ex_num : (String) Example Number
%
% % Outputs.
%
% fxy : (Matrix) Coefficient matrix of polynomial f(x,y)
%
% m : (Int) Total degree of polynomial f(x,y)



EXAMPLE_TYPE = 'From Coefficients';

switch EXAMPLE_TYPE
    case 'From Roots'
        [fxy, m] = Examples_Roots_FromRoots_Bivariate(ex_num);
        
    case 'From Coefficients'
        [fxy, m] = Examples_Roots_FromCoefficients_Bivariate(ex_num);
        
end


end



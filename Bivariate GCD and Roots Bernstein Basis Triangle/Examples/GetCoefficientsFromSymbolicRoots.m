function [fxy] = GetCoefficientsFromSymbolicRoots(root_mult_arr)
%  
%
% % Inputs.
%
% root_mult_arr : Matrix consisting of [Factor Multiplicity] pairs, where
% the factors are symbolic.
%
% % Outputs
%
% fxy : Coefficients of polynomial f(x,y)

% Get the factors of f(x,y)
arr_sym_factors_fxy = GetFactors(root_mult_arr);

% Get the coefficients of f(x,y)
fxy = GetCoefficientsFromFactors(arr_sym_factors_fxy);

end
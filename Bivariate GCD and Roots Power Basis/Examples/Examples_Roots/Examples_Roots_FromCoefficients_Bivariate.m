function [fxy, m] = Examples_Roots_FromCoefficients_Bivariate(ex_num)
%
% % Inputs
%
% ex_num : (String) Example Number
%
% % Outputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% m : (Int) Total degree of polynomial f(x,y)

% Add path to examples files
addpath(genpath('../Examples'))

% Get the symbolic roots and their multiplicities
[f_root_sym_mult_arr] = Roots_Examples_Bivariate(ex_num);

% Get the coefficients of the bivariate polynomial f(x,y)
[fxy, m, ~, ~] = GetCoefficientsFromSymbolicRoots(f_root_sym_mult_arr);





end
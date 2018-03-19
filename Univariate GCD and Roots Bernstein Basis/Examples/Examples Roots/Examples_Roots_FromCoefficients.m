function [fx, fx_root_mult_array] = Examples_Roots_FromCoefficients(ex_num)
%
% % Inputs
%
% ex_num : (String) Example Number
%
%
% % Outputs
%
% fx : (Column vector) Coefficients of the polynomial f(x) in the Bernstein form.
%
% 
%
% >> Examples_Roots_FromCoefficients('1')

% Add path to examples folder
addpath(genpath('../Examples'));

% Get the factors and corresponding multiplicities of f(x)
fx_root_mult_array = Roots_Examples_Univariate(ex_num);

% Get the coefficients of f(x) in Bernstein form
[fx] = BuildPolyFromRootsSymbolic(fx_root_mult_array);

% Get the symbolic polynomail
fx_sym = GetSymbolicPolyFromSymbolicRoots(fx_root_mult_array);

% Display the symbolic polynomial
disp(fx_sym);




end


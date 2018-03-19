function [fx, gx, dx, ux, vx] = Examples_GCD_FromCoefficients(ex_num)
%
% % Inputs
%
% ex_num : Example Number
%
% % Outputs
%
% [fx, gx] : Coefficients of the polynomials f(x) and g(x)
%
% dx : Coefficients of the common divisor d(x)
%
% [ux, vx] : Coefficients of the cofactor polynomials u(x) and v(x)

addpath(genpath('../Examples'));

% Get the rm (Root - Multiplicity) array for each of the polynomials
[f_rm_array, g_rm_array, d_rm_array, u_rm_array, v_rm_array] = ...
    GCD_Examples_Univariate_2Polys(ex_num);

% Get coefficients of polynomials f(x) and g(x)
fx = GetCoefficientsFromSymbolicRoots(f_rm_array);
gx = GetCoefficientsFromSymbolicRoots(g_rm_array);

% Get coefficients of polynomial d(x)
dx = GetCoefficientsFromSymbolicRoots(d_rm_array);

% Get coefficients of polynomials u(x) and v(x)
ux = GetCoefficientsFromSymbolicRoots(u_rm_array);
vx = GetCoefficientsFromSymbolicRoots(v_rm_array);

end



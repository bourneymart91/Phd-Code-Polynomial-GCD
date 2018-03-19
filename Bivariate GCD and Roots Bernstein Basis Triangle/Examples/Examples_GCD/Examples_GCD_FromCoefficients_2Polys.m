function [fxy, gxy, uxy, vxy, dxy, m, n, t] = Examples_GCD_FromCoefficients_2Polys(ex_num)
%
% % Inputs
%
% ex_num : Example number
%
%
% % Outputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% uxy : (Matrix) Coefficients of u(x,y)
%
% vxy : (Matrix) Coefficients of v(x,y)
%
% dxy : (Vector) Coefficients of polynomial d(x,y)
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% t : (Int) Total degree of d(x,y)

syms x y;

% Add path to examples file
addpath(genpath('../Examples'));

% Get the symbolic factors f_{i} and corresponding multiplicity m_{i}
[f_root_mult_arr,g_root_mult_arr,d_root_mult_arr,...
    u_root_mult_arr,v_root_mult_arr] = GCD_Examples_Bivariate_2Polys(ex_num);

% Get coefficients of the polynomials
[fxy] = GetCoefficientsFromSymbolicRoots(f_root_mult_arr);
[gxy] = GetCoefficientsFromSymbolicRoots(g_root_mult_arr);
[dxy] = GetCoefficientsFromSymbolicRoots(d_root_mult_arr);
[uxy] = GetCoefficientsFromSymbolicRoots(u_root_mult_arr);
[vxy] = GetCoefficientsFromSymbolicRoots(v_root_mult_arr);

% Get the symbolic polynomials in power form
symbolic_f = GetSymbolicPoly(f_root_mult_arr);
symbolic_g = GetSymbolicPoly(g_root_mult_arr);
symbolic_d = GetSymbolicPoly(d_root_mult_arr);
display(symbolic_f)
display(symbolic_g)
display(symbolic_d)


% Get the total degree of the polynomials f,g,d when in power form.
m = double(feval(symengine, 'degree', symbolic_f));
n = double(feval(symengine, 'degree', symbolic_g));
t = double(feval(symengine, 'degree', symbolic_d));

fprintf([mfilename ' : ' sprintf('Total Degree of f(x,y) : %i \n',m)]);
fprintf([mfilename ' : ' sprintf('Total Degree of g(x,y) : %i \n',n)]);
fprintf([mfilename ' : ' sprintf('Total Degree of d(x,y) : %i \n',t)]);

end




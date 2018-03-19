function [fxy, gxy, hxy, uxy, vxy, wxy, dxy,...
    m, m1, m2,...
    n, n1, n2,...
    o, o1, o2,...
    t, t1, t2] = ...
    Examples_GCD_FromCoefficients_3Polys(ex_num)
%
% % Inputs
%
% ex_num : (String) Example number
%
%
% % Outputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of polynomial h(x,y)
%
% dxy : (Matrix)Coefficients of polynomial d(x,y)
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% o : (Int) Total degree of h(x,y)
%
% t : (Int) Total degree of d(x,y)

syms x y;

addpath(genpath('../Examples'));

[f_root_mult_arr, g_root_mult_arr, h_root_mult_arr, d_root_mult_arr,...
    u_root_mult_arr, v_root_mult_arr, w_root_mult_arr] = GCD_Examples_Bivariate_3Polys(ex_num);


% Given the set of symbolic factors and corresponding multiplicities, get
% the matrices of coefficients of f(x,y), g(x,y) and h(x,y) as well as
% d(x,y), u(x,y), v(x,y) and w(x,y)
[fxy] = GetCoefficientsFromSymbolicRoots(f_root_mult_arr);
[gxy] = GetCoefficientsFromSymbolicRoots(g_root_mult_arr);
[hxy] = GetCoefficientsFromSymbolicRoots(h_root_mult_arr);

[dxy] = GetCoefficientsFromSymbolicRoots(d_root_mult_arr);

[uxy] = GetCoefficientsFromSymbolicRoots(u_root_mult_arr);
[vxy] = GetCoefficientsFromSymbolicRoots(v_root_mult_arr);
[wxy] = GetCoefficientsFromSymbolicRoots(w_root_mult_arr);

% Get the symbolic polynomials in the power basis for printing to screen.
symbolic_f = GetSymbolicPoly(f_root_mult_arr);
symbolic_g = GetSymbolicPoly(g_root_mult_arr);
symbolic_h = GetSymbolicPoly(h_root_mult_arr);

symbolic_d = GetSymbolicPoly(d_root_mult_arr);

symbolic_u = GetSymbolicPoly(u_root_mult_arr);
symbolic_v = GetSymbolicPoly(v_root_mult_arr);
symbolic_w = GetSymbolicPoly(w_root_mult_arr);


display(symbolic_f)
display(symbolic_g)
display(symbolic_h)
display(symbolic_d)


% Get the total degree of the polynomials f,g,d when in power form.
[m, m1, m2] = GetDegreeStructure(symbolic_f);
[n, n1, n2] = GetDegreeStructure(symbolic_g);
[o, o1, o2] = GetDegreeStructure(symbolic_h);
[t, t1, t2] = GetDegreeStructure(symbolic_d);

% Print Degree structure to screen
fprintf([mfilename ' : ' sprintf('Total Degree of f(x,y) : %i \n',m)]);
fprintf([mfilename ' : ' sprintf('Total Degree of g(x,y) : %i \n',n)]);
fprintf([mfilename ' : ' sprintf('Total Degree of h(x,y) : %i \n',o)]);
fprintf([mfilename ' : ' sprintf('Total Degree of d(x,y) : %i \n',t)]);

end

function [m,m1,m2] = GetDegreeStructure(symbolic_poly)
%
% % Inputs
%
% symbolic_poly : (Symbolic)
%
% % Outputs
%
% m : (Int) Degree of polynomial 
%
% m1 : (Int) Degree of symbolic polynomial with respect to x
%
% m2 : (Int) Degree of symbolic polynomial with respect to y
%

syms x y

% Get total degree
m = double(feval(symengine, 'degree', symbolic_poly));

% Get relative degree with respect to x
m1 = double(feval(symengine, 'degree', symbolic_poly,x));

% Get relative degree with respect to y
m2 = double(feval(symengine, 'degree', symbolic_poly,y));
end


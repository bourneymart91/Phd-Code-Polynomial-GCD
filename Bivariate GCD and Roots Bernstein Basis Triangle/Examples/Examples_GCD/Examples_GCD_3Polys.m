function [fxy, gxy, hxy, uxy, vxy, wxy, dxy,...
    m, m1, m2,...
    n, n1, n2,...
    o, o1, o2,...
    t, t1, t2]  = Examples_GCD_3Polys(ex_num)
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
% [fxy, gxy, hxy] : Coefficient matrix of polynomials f(x,y) g(x,y) and
% h(x,y)
%
% [uxy, vxy, wxy] : Coefficient matrix of polynomials u(x,y) v(x,y) and
% w(x,y)
%
% dxy : Coefficient matrix of polynomial d(x,y)
%
% [m, m1, m2] : Degree structure of polynomial f(x,y)
%
% [n, n1, n2] : Degree structure of polynomial g(x,y)
%
% [o, o1, o2] : Degree structure of polynomial h(x,y)
%
% [t, t1, t2] : Degree structure of polynomial d(x,y)



[fxy, gxy, hxy, ...
    uxy, vxy, wxy,...
    dxy,...
    m, m1, m2,...
    n, n1, n2,...
    o, o1, o2,...
    t, t1, t2] = Examples_GCD_FromCoefficients_3Polys(ex_num);


ex_num_var = ex_num(end);



boolBad_scaling = true;

switch boolBad_scaling
    case true
        f_multiplier = 10^5;
        g_multiplier = 10^5;
        h_multiplier = 10^(-5);
        
        LineBreakLarge();
        fprintf('BAD SCALING IS ON')
        LineBreakLarge();
        
    case false
        f_multiplier = 1;
        g_multiplier = 1;
        h_multiplier = 1;
end


switch ex_num_var
    
    case 'a'
        fxy = fxy .* f_multiplier;
        gxy = gxy .* g_multiplier;
        hxy = hxy .* h_multiplier;
        
    case 'b'
        
        fxy = fxy .* g_multiplier;
        gxy = gxy .* f_multiplier;
        hxy = hxy .* h_multiplier;
        
    case 'c'
        
        fxy = fxy .* h_multiplier;
        gxy = gxy .* f_multiplier;
        hxy = hxy .* g_multiplier;
        
end




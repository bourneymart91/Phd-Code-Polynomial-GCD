function [u,v,w,res,cond] = o_gcd_zeng(fx,gx)
% Given two polynomials, get the gcd using zeng method
% Note Zengs method takes strings as inputs, so must first convert to
% string from vector form

addpath('NAClabForMatlab2013');

 
% Get symbolic expression for f(x)
sym_fx = GetSymbolicFromCoefficients(fx);

% Get symbolic expression for g(x)
sym_gx = GetSymbolicFromCoefficients(gx);

% Get f(x) as string
string_fx = char(expand(sym_fx));

% Get g(x) as string
string_gx = char(expand(sym_gx));

 
[u,v,w,res,cond] = PolynomialGCD(string_fx,string_gx);

%[u2,v2,w2] = mvGCD(f_mat,g_mat)

end
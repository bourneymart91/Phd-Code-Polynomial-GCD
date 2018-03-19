function [f] = GetSymbolicPolyFromSymbolicRoots(sym_root_mult_arr_fxy)
% GetSymbolicPolyFromSymbolicRoots(factors_f)
%
% Given the array of factors of f(x,y) and their multiplicities, compute
% the symbolic polynomial.
%
% Inputs
%
% sym_root_mult_arr_fxy : Symbolic Factors and multiplicities of polynomial f(x,y)
%
% Outputs.
%
% f : Symbolic polynomial f(x,y)

syms x y

% Initialise symbolic polynomial
f = 1;

% Get number of distinct factors
nFactors = size(sym_root_mult_arr_fxy,1);

for i = 1:1:nFactors

    % Get the symbolic factor
    factor = sym_root_mult_arr_fxy(i,1);

    % Get the multiplicity of the factor
    mult = sym_root_mult_arr_fxy(i,2);

    % Compute f
    f = f * (factor.^mult); 
end



end
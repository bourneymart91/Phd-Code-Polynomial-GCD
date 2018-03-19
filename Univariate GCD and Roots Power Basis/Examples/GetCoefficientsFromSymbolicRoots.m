function [fx] = GetCoefficientsFromSymbolicRoots(f_root_mult_arr)
% Get the coefficients of the polynomial f(x) given an array of factors and
% their corresponding multiplicity in f(x).
%
% Inputs.
%
% f_root_mult_arr : (Matrix) Each row consists of a symbolic factor, and
% the multiplicity of the factor.
%
% Outputs.
%
% fx : Coefficients of the polynomial f(x).

syms x

% Get the number of factors
nFactors = size(f_root_mult_arr,1);

% Initialise the symbolic polynomial
sym_poly = 1;

for i = 1:1:nFactors
   
    % Get the factor
    sym_factor = f_root_mult_arr(i,1);
    
    % Get the multiplicity of the factor in f(x)
    mult = f_root_mult_arr(i,2); 
    
    sym_factor = sym_factor^mult;
        
    sym_poly = sym_poly * sym_factor;
end

% Get coefficients of the polynomial 
fx = double(fliplr(coeffs(sym_poly,x,'All')))';
end
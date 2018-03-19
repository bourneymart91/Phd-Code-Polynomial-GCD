function [fxy, m, m1, m2] = GetCoefficientsFromSymbolicRoots(f_root_mult_arr)
%
% % Inputs
%
% f_root_mult_arr : (Matrix) Each row contains a symbolic factor and the
% multiplicity of the factor in f(x,y)
%
% % Outputs%
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) 
%
% m : (Int) Total degree of f(x,y)
%
% m1 : (Int) Degree of f(x,y) with respect to x
%
% m2 : (Int) Degree of f(x,y) with respect to y

% Start symbolic characters
syms x y

% Get number of distinct factors in f(x,y)
nFactors = size(f_root_mult_arr,1);


% Initialise the symbolic polynomial
temp_sym_poly = 1;

for i = 1:1:nFactors
   
    % Get the symbolic factor
    sym_factor = f_root_mult_arr(i,1);
    
    % Get the multiplicity of the factor
    mult = f_root_mult_arr(i,2); 
    
    % Get expression for factor expanded
    sym_factor = sym_factor^mult;
        
    % Multiply by temp polynomial.
    temp_sym_poly = temp_sym_poly * sym_factor;
end


symbolic_poly = temp_sym_poly;
display(symbolic_poly)

% Get coefficients of the polynomial 
fxy = double(rot90(coeffs(symbolic_poly,[x,y],'All'),2));

% Get degree structure of f(x,y)
m = double(feval(symengine, 'degree', symbolic_poly));
m1 = double(feval(symengine, 'degree', symbolic_poly,x));
m2 = double(feval(symengine, 'degree', symbolic_poly,y));

end
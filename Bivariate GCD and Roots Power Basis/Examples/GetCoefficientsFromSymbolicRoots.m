function [fxy,m,m1,m2] = GetCoefficientsFromSymbolicRoots(f_root_mult_arr)
%
%
%

syms x y

nFactors = size(f_root_mult_arr,1);


% Initialise the symbolic polynomial
temp_sym_poly = 1;

for i = 1:1:nFactors
   
    % Get the factor
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

m = double(feval(symengine, 'degree', symbolic_poly));
m1 = double(feval(symengine, 'degree', symbolic_poly,x));
m2 = double(feval(symengine, 'degree', symbolic_poly,y));

end

function fx = GetSymbolicPolyFromSymbolicRoots(f_factor_mult_array)

% This function takes the factors of the polynomial f(x) in a symbolic
% form, gets the coefficients of the factors in Bernstein form, then
% multiplies these to form polynomial f(x) in bernstein form.

syms x

% Get the number of unique factors
nFactors = size(f_factor_mult_array,1);

% Initialise polynomial f(x)
fx = 1;

for i = 1 : 1 : nFactors
   
   % Get the factor
   sym_factor = f_factor_mult_array(i,1);
   
   % Get the multiplicity
   mult = f_factor_mult_array(i,2);
   
   fx = fx .* (sym_factor.^(mult));
   
end



end
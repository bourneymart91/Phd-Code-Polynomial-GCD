function [fx] = BuildPolyFromRootsSymbolic(f_factor_mult_array)
% This function takes the factors of the polynomial f(x) in a symbolic
% form, gets the coefficients of the factors in Bernstein form, then
% multiplies these to form polynomial f(x) in bernstein form.
%
% % Inputs
%
% f_factor_mult_array : (Matrix) Where first column consists of factors and
% second column consists of their multiplicity in f(x)
%
% % Outputs
%
% fx : (Vector) Coefficients of polynomial f(x)


syms x

nFactors = size(f_factor_mult_array, 1);

for i = 1 : 1 : nFactors
    
    % Get the factor
    sym_factor = f_factor_mult_array(i,1);
    
    % Get the multiplicity
    mult = f_factor_mult_array(i,2);
    
    % Get the expression (factor).^{mult}
    sym_factor = sym_factor^mult;
    
    % Get coefficients in power basis
    try
        pwr_poly = double(fliplr(coeffs(sym_factor,x,'All')))';
    catch
        pwr_poly = double((coeffs(sym_factor,x)))';
    end
    
    % Convert the factor to a polynomial in Bernstein form
    arr_factors_fx{i} = PowerToBernstein(pwr_poly);
    
end

% Multiply all matrices of coefficients of factors to form coefficients of
% polynomial f(x)

fx = arr_factors_fx{1};

for i = 2:1:nFactors
    fx = Bernstein_Multiply(fx,arr_factors_fx{i});
    
end

end
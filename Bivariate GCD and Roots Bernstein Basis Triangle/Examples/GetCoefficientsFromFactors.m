function [fxy] = GetCoefficientsFromFactors(arr_sym_factors_fxy)
% GetCoefficientsFromFactors(arr_sym_factors_fxy)
%
% Given an array of factors of polynomial f(x,y) compute the coefficients
% of f(x,y) in Bernstein form.
%
% % Inputs.
%
% arr_sym_factors_fxy : 
%
% % Outputs
% 
% fxy : Coefficients of polynomial f(x,y)


syms x y;


nFactors = length(arr_sym_factors_fxy);
arr_factors = cell(nFactors,1);

% For each factor
for i = 1:1:length(arr_sym_factors_fxy)
    
    % Get the factor in power form
    factor = double(rot90(coeffs(arr_sym_factors_fxy{i},[x,y],'All'),2));
    
    % Get total degree of the factor
    arr_m{i} = double(feval(symengine, 'degree', arr_sym_factors_fxy{i}));
    
    % Get the factor in the Bernstein form
    factor_brn = PowerToBernstein(factor, arr_m{i});
    arr_factors{i} = factor_brn;
    
    
end

% Get product of all factors
fxy = arr_factors{1};
m = arr_m{1};

for i = 2:1:length(arr_factors)
    
    % m = degree of f(x)
    
    % n = degree of next factor
    n = arr_m{i};
    
    fxy = Bernstein_Multiply(fxy,arr_factors{i},m,n);
    
    m = m + n;
    
end



end
function [fx] = Examples_Roots(ex_num)


EXAMPLE_TYPE = 'From Coefficients';

switch EXAMPLE_TYPE
    case 'From Roots'
        
        % Get a set of roots and multiplicities for f(x)
        fx_root_mult_array = Examples_Roots_FromRoots(ex_num);
        
        % Get the coefficients of the polynomial f(x)
        fx = GetCoefficientsFromRoots(fx_root_mult_array);
        
        % Print the roots and coefficients of f(x)
        PrintFactorization(fx_root_mult_array,'f');
        
        PrintCoefficientsBivariate(fx,'f');
        
        
        
    case 'From Coefficients'
        
        % Get the coefficients of polynomial f(x)
        fx = Examples_Roots_FromCoefficients(ex_num);
        PrintCoefficientsBivariate(fx,'f')
        
    otherwise
        error('Example polynomial is either from roots or from coefficients')
end

end


function [factor_mult_arr] = Deconvolution_Examples_Univariate(ex_num)
% Get the factors and multiplicity of each factor in f(x,y)
%
% % Inputs
%
% ex_num : Example Number
%
% % Outputs
%
% factor_mult_arr : Symbolic factors and corresponding multiplicity

% Input f_{i} polynomials
x = sym('x');

% Set example number

switch ex_num
    case '1'
        
        % Create set of factors whose multiplicities are defined in vMult
        factor_mult_arr = ...
            [
                (x + 0.017746571505)  20;
                (x - 0.529678501354)  40;
            ];
        
        
    case '2'
        
        factor_mult_arr = ...
            [
                (x - 0.1259)  1
                (x - 0.2789)  3
                (x - 0.389)   4
                (x - 0.4213)  4
                (x - 0.5432)  5
                (x + 0.7923)  12
            ];
        
        
        
    case '3'
        
        factor_mult_arr = ...
            [
                (x - 2)           2
                (x - 3.2789)      4
                (x - 1.589)      12
            ];
        
        
    case '4'
        
        factor_mult_arr = ...
            [
                (x - 0.56897)     3
                (x + 1.24672)     6
                (x + 0.56921)     9
            ];
        
        
    case '5'
        factor_mult_arr = ...
            [
                (x - 0.246512)  2
                (x - 1.214654)  5
                (x + 0.567890)  7
                (x + 0.214654)  12
            ];
    case '6'
        factor_mult_arr = ...
            [
                (x - 0.246512)   2
                (x - 1.214654)   5
                (x + 0.567890)   7
                (x + 0.214654)   12
            ];
    case '7'
        
        factor_mult_arr = ...
            [
                (x-2)       1
                (x-3.2789)  3
                (x-1.589)   4
                (x-0.7213)  4
                (x-1.5432)  5
                (x+5.72)    12
            ];
        
    case '8'
        
        factor_mult_arr = ...
            [
                (x - 2.1234565487)       1
                (x - 1.589212457)   4
                (x - 0.7213)  10
                (x - 6.5432)  7
                (x + 0.72)    20
            ];
        
    case '9'
        
        factor_mult_arr = ...
            [
                (x - 2.16547697898)         2
                (x + 0.2789)    30
                (x - 1.589)     12
            ];
        
    case '10'
        
        factor_mult_arr = ...
            [
                (x - 3.654132475632154)         4
                (x + 0.278912456789)    40
                (x - 1.589)     12
            ];
        
    case '11'

    factor_mult_arr = ...
        [
            (x - 3.654132475632154)        10
            (x + 0.278912456789)    40
            (x - 1.589)     12
        ];
    
    case '12'

    factor_mult_arr = ...
        [
            (x - 0.654132475632154)        10
            (x + 0.278912456789)    40
            (x - 1.589)     4
            (x + 0.79877436257878) 5
        ];
    case '13'

    factor_mult_arr = ...
        [
            (x - 0.654132475632154)         9
            (x + 0.278912456789)            16
            (x - 1.589)                     7
            (x + 0.79877436257878)          20
        ];
            
    otherwise
        error([mfilename ' : Not a valid example number'])
end
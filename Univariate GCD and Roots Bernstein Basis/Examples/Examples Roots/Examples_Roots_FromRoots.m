function [root_mult_array_fx,m] = Examples_Roots_FromRoots(ex_num)
% Examples_Roots(n)
%
% Get a set of roots and multiplicities of a polynomial f(x).
%
% Inputs.
%
%
% n : Example Number
%
% Outputs
%
% root_mult_array_fx : Matrix of roots of f(x) and their corresponding
%                      multiplicities. Each row is a root, multiplicity
%                      pair.
%
% m : Degree of polynomial f(x).

pattern = 'Custom:m=(\d+).low=(-?\d+).high=(-?\d+)';

if ~isempty(regexp(ex_num,pattern,'start'))
  
    str = ex_num;
    
    expression_m = regexp(str,'m=(\d+)','tokens');
    m_str = expression_m{1};
     
    expression_low = regexp(str, 'low=(-?\d+)', 'tokens');
    low_str = expression_low{1};
    
    expression_high = regexp(str, 'high=(-?\d+)','tokens');
    high_str = expression_high{1};
    

    
    m = str2double(m_str);    
    intvl_low = str2double(low_str);
    intvl_high = str2double(high_str);
            
    root_mult_array_fx = BuildRandomPolynomial(m,intvl_low,intvl_high);
    display(root_mult_array_fx);
    
    
else
    
    
    
    
    switch ex_num
        
        case 'Example'
            root_mult_array_fx =...
                [
                0.1     7;
                0.9     12;
                ];
            
        case 'Example Zeng'
            root_mult_array_fx = ...
                [
                10/11   5
                20/11   3
                30/11   2
                ];
            
        case '-1'
            root_mult_array_fx = ...
                [
                0.1     1;
                0.5     10;
                ];
            
        case '0'
            root_mult_array_fx = ...
                [
                0.1     1;
                0.7     1;
                0.9     1;
                ];
            
        case '1'
            root_mult_array_fx = ...
                [
                0.1     1;
                0.9     2;
                0.2     3;
                0.5     4;
                1.6     5;
                ];
            
        case '2'
            root_mult_array_fx = ...
                [
                0.5     1;
                0.8     2;
                ];
            
        case '3'
            root_mult_array_fx = ...
                [
                0.3,     1;
                0.6,    2;
                -1.1,    3;
                ];
        case '4'
            root_mult_array_fx = ...
                [
                0.1    1;
                0.5    2;
                0.9    3;
                1.4    4;
                ];
            
        case '5'
            root_mult_array_fx = ...
                [
                0.1   1;
                0.5   2;
                0.8   3;
                0.4   4;
                1.3   5;
                ];
            
        case '6'
            root_mult_array_fx = ...
                [
                0.14,  1;
                0.56,  2;
                0.89,  3;
                0.45,  4;
                0.37,  5;
                1.3,  6
                ];
            
        case '7'
            root_mult_array_fx = ...
                [
                0.14    1;
                0.56    2;
                0.89    3;
                1.45    4;
                2.37    5;
                -3.6    6;
                0.8     7;
                ];
            
        case '8'
            root_mult_array_fx = ...
                [
                0.14,  1;
                0.56,  2;
                0.89,  3;
                1.45,  4;
                2.37,  5;
                -3.61,  6
                0.8     7
                0.6     8
                ];
            
        case '9'
            root_mult_array_fx= ...
                [
                0.14,  1;
                0.56,  2;
                0.89,  3;
                1.45,  4;
                2.37,  5;
                -3.61,  6
                0.8     7
                0.6     8
                1.2     9
                ];
        case '10'
            root_mult_array_fx = ...
                [
                0.14,  1;
                0.56,  2;
                0.89,  3;
                1.45,  4;
                2.37,  5;
                -3.61,  6
                0.8     7
                0.6     8
                1.2     9
                0.10003245657   10
                ];
            
            % cases with roots of same multiplicites
        case '11'
            root_mult_array_fx = ...
                [
                0.1     1;
                0.4     2;
                0.7     2;
                ];
            
        case '12'
            root_mult_array_fx = ...
                [
                0.1     1;
                0.5     2;
                0.6     2;
                0.9     3;
                1.1     3;
                ];
            
        case '13'
            root_mult_array_fx = ...
                [
                0.1     1;
                0.5     2;
                0.6     2;
                0.9     3;
                1.1     5;
                ];
            
        case '14'
            root_mult_array_fx = ...
                [
                0.1     1;
                0.5     2;
                0.6     2;
                0.9     3;
                1.1     7;
                ];
            
            
        case '15'
            root_mult_array_fx = ...
                [
                3.654132475632154       10
                -0.278912456789         40
                1.589                   12
            ];
            
        case 'Custom'
            
            prompt = 'Enter the degree of Polynomial f(x) :';
            m = input(prompt);
            
            intvl_low = -1;
            intvl_high = +1;
            
            root_mult_array_fx = BuildRandomPolynomial(m,intvl_low,intvl_high);
        otherwise
            error('Example Number not found')
            
    end
    
end

m = sum(root_mult_array_fx(:,2)-1);

end
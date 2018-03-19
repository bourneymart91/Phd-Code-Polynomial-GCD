function root_mult_array_fx = Examples_Roots_FromRoots(ex_num)
% Given an example number, return a polynomial in terms of its roots and
% their corresponding multiplicities.
%
% Input
%
% ex_num : Example Number
%
% Outputs.
%
% root_mult_array_fx : Array of roots and corresponding multiplicities in f(x)



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
            root_mult_array_fx = ...
                [
                1   1
                2   2
                3   3
                5   5
                ];
            
        case '1'
            root_mult_array_fx = ...
                [
                0.1 7;
                0.9 12;
                ];
        case '2'
            root_mult_array_fx = ...
                [
                0.1 1;
                0.2 2;
                0.9 3;
                ];
            
        case '3'
            root_mult_array_fx = ...
                [
                6   10;
                -3  3;
                -7  1;
                ];
        case '4'
            root_mult_array_fx = ...
                [
                1   1
                2   2
                3   3
                4.5 1
                ];
        case '5'
            root_mult_array_fx = ...
                [
                1.05467 1
                2.24587 2
                5.54743 3
                1.75647 2
                ];
        case '6'
            root_mult_array_fx = ...
                [
                1.05467 1
                2.24587 2
                5.54743 3
                1.75647 2
                2.56478 5
                1.15445 1
                ];
        case 'Custom'
            prompt = 'Enter the degree of Polynomial f(x) :';
            m = input(prompt);
            
            intvl_low = -1;
            intvl_high = +1;
            
            root_mult_array_fx = BuildRandomPolynomial(m,intvl_low,intvl_high);
            
        case 'Wilkinsons'
            root_mult_array_fx = ...
                [
                1 1
                2 1
                3 1
                4 1
                5 1
                6 1
                ];
        case 'Zeng'
            % % A Method Computing Multiple Roots of Inexact Polynomials
            root_mult_array_fx = ...
                [
                1   20
                2   15
                3   10
                4   5
                ];
            
        case 'Zeng2'
            % A Method Computing Multiple Roots of Inexact Polynomials - Zeng
           
            root_mult_array_fx = ...
                [
                    1   4
                    2   3
                    3   2
                    4   1
                ];
        case 'Zeng3'
            % A Method Computing Multiple Roots of Inexact Polynomials - Zeng
            root_mult_array_fx = ...
                [
                    10/11   5
                    20/11   3
                    30/11   2
                    ];
                
        case 'Zeng4'
            % A Method Computing Multiple Roots of Inexact Polynomials - Zeng
            e = 0.1;
            root_mult_array_fx = ...
                [
                    1-e     20;
                    1       20;
                    -0.5    5;
                ];
            
        otherwise
            error('err')
    end
    
    
end

end


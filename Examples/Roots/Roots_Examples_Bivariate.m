function [factor_mult_arr] = Roots_Examples_Bivariate(ex_num)
%
% % Inputs
%
% ex_num : Example Number as a string
%
% % Outputs
%
% factor_mult_arr : Matrix where each row contains a symbolic factor and
% its multiplicity in f(x,y).

syms x y;

switch ex_num
    
    case '1'
        
        % Polynomial sequence
        % f_{0} : t = 11 , t1 = 11, t2 = 11
        % f_{1} : t = 9  , t1 = 9,  t2 = 9
        % f_{2} : t = 7  , t1 = 7,  t2 = 7
        % f_{3} : t = 5  , t1 = 5,  t2 = 5
        % f_{4} : t = 3  , t1 = 3,  t2 = 3
        % f_{5} : t = 2  , t1 = 2,  t2 = 2
        % f_{6} : t = 1  , t1 = 1,  t2 = 1
        % f_{5} : t = 0  , t1 = 0,  t2 = 0
        
        % f_{0} = (x+y-0.5)^4 * (x+y-0.75)^7
        % f_{1} = (x+y-0.5)^3 * (x+y-0.75)^6
        % f_{2} = (x+y-0.5)^2 * (x+y-0.75)^5
        % f_{3} = (x+y-0.5)^1 * (x+y-0.75)^4
        % f_{4} =               (x+y-0.75)^3
        % f_{5} =               (x+y-0.75)^2
        % f_{6} =               (x+y-0.75)^1
        
        
        factor_mult_arr = [...
            (x + y - 0.5)       4 
            (x + y - 0.75)      7
            ];
    
    case '2'
        % Polynomial sequence
        % f_{0} : t = 14 , t1 = 7, t2 = 7
        % f_{1} : t = 11 , t1 = 5, t2 = 6
        % f_{2} : t = 8  , t1 = 3, t2 = 5
        % f_{3} : t = 5  , t1 = 1, t2 = 4 
        % f_{4} : t = 3  , t1 = 0, t2 = 3
        
        
        factor_mult_arr = [...
            (x - 1.5)   3 
            (y + 0.75)  7 
            (x - 10.1)  3
            ];
        
    case '3'
        
        factor_mult_arr = [...
            (x - 1.5)       4 
            (x + y + 0.75)      7
            (x - 10.1)      3
            (x + y -0.17523547)  5
            ];
    case '4'
        factor_mult_arr = [...
            (x^2 + y^2 + 0.5) 10
            (x + 1) 3
        ];      
            
    case '5'
            
      factor_mult_arr = [...
            (x + y - 0.5)       6
            (x + y - 0.75)      9
            (x - 0.457)         3
            ];
        
        
    case '6'
        
        factor_mult_arr = [...
            (x + y - 0.5)       7
            (x + y - 0.75)      5
            (x - 0.457)         10
            (y - 0.14567987)    5
            ];
        
        
    otherwise 
        error('err')

end


end
function CP = Examples_3DBezier_ControlPoints(ex_num)
% Given an example number, return the control points of a planar Bezier 
% curve.

switch ex_num
    
    case '1'
        CP = ...
            [
                1   5   10  12      ;
                1   7   9   4       ;
                1   8   2   6       ;
            ];
    case '2'
        CP = ...
        [
             -2   0   1 ;        
             -1   0   2 ;       
              1  -5   7 ; 
        ];
    otherwise 
        error('err')
        

end
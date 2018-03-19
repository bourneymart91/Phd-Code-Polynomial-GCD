function [C1,C2] = Examples_Explicit_Bezier_Intersection(ex_num)


switch ex_num
    case '1'
        C1 = [...
            0   0.25    0.5     0.75    1 ;
            0   2       1       3       -1 ;
        ];

        C2 = [...
            0       0.25    0.5     0.75    1 ;
            4.2     -3       -3      5      4.2;
        ];
    case '2'
        
        C1 = [...
            0   0.25    0.5     0.75    1 ;
            0   2       1       3       -1 ;
        ];
    
        C2 = [...
            0   0.25    0.5     0.75    1 ;
            1.575   1.575       1.575       1.575       1.575 ;
        ];
end

end
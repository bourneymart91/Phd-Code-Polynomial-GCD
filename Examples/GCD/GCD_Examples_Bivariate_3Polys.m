function [f, g, h, d_rm, u_rm, v_rm, w_rm] = GCD_Examples_Bivariate_3Polys(ex_num)
%
% % Inputs
%
% ex_num : (String) Example number and variant
%
% % Outputs
%
% f : (Matrix)
%
% g : (Matrix)
%
% h : (Matrix)
%
% d : (Matrix)
%
% u : (Matrix)
%
% v : (Matrix)
%
% w : (Matrix)

syms x y;



ex_num_variant = ex_num(end);
ex_num_number = ex_num(1:end-1);

[d_rm, u_rm, v_rm, w_rm] = ...
    GetPolys(ex_num_number, ex_num_variant);

f = [u_rm; d_rm];
g = [v_rm; d_rm];
h = [w_rm ; d_rm];

end


function [d_rm, u_rm, v_rm, w_rm, d_root_mult_arr] = GetPolys(ex_num, ordering)


syms x y


common_1_2_rm = [];
common_1_3_rm = [];
common_2_3_rm = [];

switch ex_num
    
    case '1'
        d_root_mult_arr = [...
            (x - 0.72)          1
            (y - 0.15)          1
            (y^2 - 1.7)         1
            ];
        
        Poly1_rm = [...
            (x + 0.75)          1
            (y - 0.75)          2
            ];
        
        Poly2_rm= [...
            (x^2 + y^2 + 0.7)   1
            (x - 0.192)         1
            ];
        
        Poly3_rm= [...
            x^2 + y^2 - 0.34    3
            x - 1.91987         4
            ];
        
        
        common_1_2_rm = [...
            (x - 0.52)          2
            (x + y - 0.5 )      5
            ];
        
    case '2' % Example 8.3 in my report
        
        Poly1_rm = [...
            (x + y + 31/2500)   5
            ];
        
        Poly2_rm = [...
            (x + y + 282/625)   3
            ];
        
        Poly3_rm = [...
            12*x^2 + y^2 - 52.34    1
            ];
        
        d_root_mult_arr = [...
            (x + 14/25)             1
            (x^2 + y^2 + 51/100)    2
            (x + y + 28/25)         3
            ];
        
        
        
    case '3'
        
        Poly1_rm = ...
            [
            (y - 1/5)     2
            ];
        
        Poly2_rm = ...
            [
            (x - 3/10)    2
            (y - 2/5)     1
            ];
        
        
        Poly3_rm = [...
            12*x^2 + y^2 - 52.34    1
            ];
        
        d_root_mult_arr = ...
            [
            (x+1/10)    4
            (x-2/5)     10
            ];
        
        
        
        
    case '4'
        
        Poly1_rm = [...
            (y-0.2)     2
            ];
        
        Poly2_rm = [...
            (y-0.4)     1
            (x - 0.3)   3
            ];
        
        Poly3_rm = [...
            12*x^2 + y^2 - 52.34    1
            ];
        
        d_root_mult_arr = [...
            (x + 1)     1
            (x + 0.8)   4
            ];
        
        
        
        
    case '5'
        
        Poly1_rm = [
            (x + y + 3.0124)    1
            ];
        
        Poly2_rm = [
            (x + y + 5.4512)    1
            ];
        
        Poly3_rm = [...
            12*x^2 + y^2 - 52.34    1
            ];
        
        d_root_mult_arr = [
            (x + y + 1)         1
            (x + y + 2)         1
            ];
        
        
    case '6'
        
        Poly1_rm = [
            (x + y + 3.0124)    2
            ];
        
        Poly2_rm = [
            (x + y + 5.4512)    1
            ];
        
        Poly3_rm = [...
            12*x^2 + y^2 - 52.34    1
            ];
        
        d_root_mult_arr = [
            (x + y + 1)         1
            (x + y + 2)         2
            (x+1)               3
            ];
        
        
    case '7'
        
        Poly1_rm = [
            (x + y - 0.0124)    6
            (y - 0.789561573)   2
            ];
        
        Poly2_rm = [
            (x + y - 0.4512)    3
            (x - 1.1654987)     2
            ];
        
        Poly3_rm = [...
            12*x^2 + y^2 - 0.346587261657    1
            ];
        
        d_root_mult_arr = [
            (x^2 + y^2 - 0.51654654)    2
            (x - y - 1.124548798)       3
            (x - 0.56)                  1
            ];
        
    case '8'
        
        d_root_mult_arr = [
            (x + y + 0.0124)    6
            (x + y + 0.923547)  2
            (x + y - 0.123456)  4
            ];
        Poly1_rm = [
            (x-1)               1
            (x-2)               2
            ];
        Poly2_rm = [
            (x^2 - 5*x)         1
            ];
        
        Poly3_rm = [...
            12*x^2 + y^2 - 52.34    1
            ];
        
        
    case '9'
        
        d_root_mult_arr = [
            (x-0.5)     1
            (x-0.2)     2
            (x-0.3)     3
            (y-0.5)     6
            ];
        
        Poly1_rm = [
            (x-0.4445)  4
            (y-0.4)     4
            ];
        
        Poly2_rm = [
            (x-0.1)     5
            (y-0.2234)  5
            ];
        
        Poly3_rm = [...
            12*x^2 + y^2 - 52.34    1
            ];
        
        
    case '10'
        d_root_mult_arr = [
            (x-0.5)     1
            (x-0.2)     2
            (x-0.3)     3
            ];
        
        Poly1_rm = [
            (x-0.4445)  4
            ];
        
        Poly2_rm = [
            (x-0.1)    5
            ];
        
        Poly3_rm = [...
            12*x^2 + y^2 - 52.34    1
            ];
        
        
        
    case '11'
        d_root_mult_arr = [
            (x + 14/25) 1
            (x^2 + y^2 + 51/100)    2
            (x + y + 28/25) 3
            ];
        
        Poly1_rm = [
            (x + y + 31/2500) 6
            ];
        
        Poly2_rm = [
            (x + y + 282/625) 3
            ];
        
        Poly3_rm = [...
            12*x^2 + y^2 - 52.34    1
            ];
        
        
        
    case '12'
        
        d_root_mult_arr = [
            (x + 2.21657951321) 1
            (x^2 + y^2 + 0.5679814324687)    2
            (x + y + 42.46578784351654) 3
            ];
        
        Poly1_rm = [
            (x + y + 31/2500) 6
            (x - 0.554687987932164654)  3
            ];
        
        Poly2_rm = [
            (x + y + 282/625) 3
            ];
        
        Poly3_rm = [...
            12*x^2 + y^2 - 52.34    1
            ];
        
        
        
    case '13'
        
        d_root_mult_arr = [
            (x + 0.21657951321)                 1
            (x^2 + y^2 + 0.5679814324687)       2
            (x + y - 0.46578784351654)          3
            ];
        
        Poly1_rm = [
            (x + y + 0.0124)                    6
            (x - 0.554687987932164654)          3
            ];
        
        Poly2_rm = [
            (x + y + 0.4512)                6
            (x^2 + y^2 - 0.00104751807)     3
            ];
        
        Poly3_rm = [...
            12*x^2 + y^2 - 0.348798             1
            ( y - 0.2465879841351465498)        4
            ];
        
        
    case '14'
        
        d_root_mult_arr = [
            (x + 0.21657951321)                 1
            (x^2 + y^2 + 0.5679814324687)       2
            (x + y - 0.46578784351654)          3
            (x^2 + y^2 + 1.657981435498798)       2
            ];
        
        Poly1_rm = [
            (x + y + 31/2500)               6
            (x + y^2 + 31/2500)               6
            (x - 0.554687987932164654)      3
            ];
        
        Poly2_rm = [
            (x + y + 282/625)               6
            (x^2 + y^2 - 0.00104751807)  3
            ];
        
        Poly3_rm = [...
            (12*x^2 + y^2 - 0.348798)             1
            ( y - 0.2465879841351465498)        4
            ];
        
        
    case '15'
        
        
        Poly1_rm = [
            (x - 0.7576)                1
            (x + y - 0.0124)            1
            (x + 2*y - 0.6)             2
            (x + 1.2*y - 0.6578654)     5
            ];
        
        Poly2_rm = [
            (x + y - 0.4512)        1
            (x - 0.4512)            3
            ];
        
        Poly3_rm = [...
            0.12*x^2 + y^2 - 2.345456   1
            0.12*x^2 + y^2 - 0.34654    4
            (x - 0.701657984)           1
            ];
        
        d_root_mult_arr = [
            (x + y - 0.1654987)         4
            (x + y - 0.7454987)         3
            (x + y^2 -0.9873654)        2
            (x + 0.5)                   3
            ];
        
        
    otherwise
        error([mfilename ' : error : Not a valid example number'])
end



switch ordering
    
    case 'a'
        d_rm = d_root_mult_arr;
        u_rm = [...
            Poly1_rm; ...
            common_1_2_rm; ...
            common_1_3_rm...
            ];
        
        v_rm = [...
            Poly2_rm; ...
            common_1_2_rm; ...
            common_2_3_rm ...
            ];
            
        w_rm = [...
            Poly3_rm; ...
            common_1_3_rm; ...
            common_2_3_rm ...
            ];
        
    case 'b'
        
        d_rm = d_root_mult_arr;
        u_rm = [...
            Poly2_rm; ...
            common_1_2_rm; ...
            common_2_3_rm ...
            ];
        v_rm = [...
            Poly1_rm; ...
            common_1_2_rm; ...
            common_1_3_rm ...
            ];
        w_rm = [...
            Poly3_rm; ...
            common_1_3_rm; ...
            common_2_3_rm ...
            ];
        
    case 'c'
        
        d_rm = d_root_mult_arr;
        u_rm = [...
            Poly3_rm; ...
            common_1_3_rm; ...
            common_2_3_rm ...
            ];
        v_rm = [...
            Poly2_rm; ...
            common_1_2_rm; ...
            common_2_3_rm ...
            ];
        w_rm = [...
            Poly1_rm; ...
            common_1_2_rm; ...
            common_1_3_rm...
            ];
    otherwise
        error('Not valid ordering')
end


end
function [fxy_matrix, gxy_matrix] = Examples_Intersections(ex_num)


f_roots_x = [];
f_roots_y = [];
f_roots_xy = [];
g_roots_x = [];
g_roots_y = [];
g_roots_xy = [];

switch ex_num
    case '1'
        f_roots_x =...
            [
                1   1;
                2   2;
                1.5 3;
            ];
    
        g_roots_x = ...
            [
                1   1;
                2   2;
                3   3;
            ];
        
end



f_roots_x = mult_roots_x(f_roots_x);
f_roots_y = mult_roots_y(f_roots_y);
g_roots_x = mult_roots_x(g_roots_x);
g_roots_y = mult_roots_y(g_roots_y);

f_roots = [f_roots_x; f_roots_y ; f_roots_xy];
g_roots = [g_roots_x; g_roots_y ; g_roots_xy];


fxy_matrix = BuildPoly_NonSeparable(f_roots)
gxy_matrix = BuildPoly_NonSeparable(g_roots)

[r,c] = size(fxy_matrix);
m1 = r -1;
m2 = c -1;

[r,c] = size(gxy_matrix);
n1 = r -1;
n2 = c -1;


end
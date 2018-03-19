
function [f_x_roots,f_y_roots] = Examples_Bivariate_Separable(ex_num)
switch ex_num
    case '1'
        f_x_roots = [
            0.8   1;
            0.5   1;
            ];
        f_y_roots = [
            0.8   1;
            0.2  1;
            ];
    case '2'
        f_x_roots = [
            0.8   1;
            ];
        f_y_roots = [
            0.2  1;
            ];
end
end

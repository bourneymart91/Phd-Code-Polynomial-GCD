function [] = o1_bivariate_Pwr()

example_style = 2;

switch example_style
    case 1  
        ex_num = 1;
        
        % Get the roots of polynomial f
        [f_x_roots,f_y_roots] = Examples_Bivariate_Separable(ex_num);
        
        % Get the coefficients of the polynomials f(x) and f(y) in the power basis.
        f_x_poly_pwr = BuildPoly_Pwr(f_x_roots)
        f_y_poly_pwr = BuildPoly_Pwr(f_y_roots)
        
        % Get the degrees m1 and m2 of polynomials f wrt x and f wrt y,
        % respectively.
        m1 = length(f_x_poly_pwr)-1;
        m2 = length(f_y_poly_pwr)-1;
        
        
        % Get the matrix of coefficients of the polynomial surface in power form.
        fxy_matrix_Pwr = f_x_poly_pwr * f_y_poly_pwr'
        fxy_matrix_Pwr
    case 2
        ex_num = 6;
        switch ex_num
            case 1
        fxy_matrix_Pwr = [  1 1 1; 
                            2 3 5;
                            8 -1 4]
            case 2
                fxy_matrix_Pwr = ...
                    [ 0   0    0      1    -5;
                      0   0    1     -7     10;
                      0   1   -6      2     15;
                      1  -7    6      28   -40;
                     -1   6   -1     -24    20
                    ];
                
            case 3
                %   BOUNDS FOR THIS PROBLEM
                %   x = [0, .. , 3.2]
                %   y = [-3, .. , 1]
                fxy_matrix_Pwr = ...
                    [  1   3  0  -4;
                      -5 -15  0  20;
                       7  21  0  -28;
                      -3  -9  0  12
                     ];
            case 4
                %
                %
                %
                fxy_matrix_Pwr = ...
                    [  1   3  0  -4;
                      -1  -3  0   4;
                     ];
            case 5
                %
                %
                %
                fxy_matrix_Pwr = ...
                    [
                    1 3 0 4
                    ];
            case 6
                % Bivariate - Non Separable - Easier Example
                fxy_matrix_Pwr = ...
                    [
                        0   1   -5;
                        1   -7  10;
                        -1  6   -5;
                    ];
        end
                        
        m1 = size(fxy_matrix_Pwr,1) -1;
        m2 = size(fxy_matrix_Pwr,2) -1;
   
        
end

plot_fxy(fxy_matrix_Pwr)
end
function [fxy, m] = Examples_Roots_FromRoots_Bivariate(ex_num)
%
% % Inputs
%
% ex_num : (String) Example Number
%
% % Outputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% m : (Int) Total degree of f(x,y)

root_mult_array_f_x = [];
root_mult_array_f_y = [];
polys_xy = [];

switch ex_num
    
    
    
    
    case '0'
        root_mult_array_f_y = ...
            [
            1   1
            2   2
            ];
        
        m = 3;
    case '1'
        % Relates to roots_examples tex file example 1
        
        % Get roots of f with respect to x
        root_mult_array_f_x =...
            [
            0.9   2;
            ];
        
        % Get roots of f with respect to y
        root_mult_array_f_y = ...
            [
            0.3   1
            ];
        
        % Get roots of f with respect to x and y
        % (x+y-0.1)
        polys_xy{1,1} = ...
            [
            -0.1  1;
            1  0
            ];
        
        m = 4;    
    case '2'
        % Example 1
        % in File Bivariate Root Finding - By example.tex
        %
        % (x-1)^2
        % (x-3)
        % (y+2)^2
        % (y-1)
        
        root_mult_array_f_x =[...
            1  2;
            3  1;
            ];
        
        root_mult_array_f_y =[...
            -2  2;
            1   1;
            ];
        
        m = 6;
        
        
    case '3'
        %         ( x - 0.7526)^{2}
        %         ( x - 0.1236)^{3}
        %         ( x - 1.1000)^{4}
        %         ( y - 0.1459)^{2}
        %         ( y - 0.9876)^{1}
        %         (-0.75 x^{0}y^{0} + 1.00 x^{1}y^{0} +1.00 x^{0}y^{1} +0.00 x^{1}y^{1} )^{1}
        %         ( 2.62 x^{0}y^{0} + 1.00 x^{1}y^{0} +1.00 x^{0}y^{1} +0.00 x^{1}y^{1} )^{2}
        
        
        polys_xy{1,1} = [...
            -0.753  1;
            1       0;
            ];
        polys_xy{2,1} = [...
            2.6235  1;
            1       0;
            ];
        polys_xy{3,1} = [...
            2.6235  1;
            1       0;
            ];
        
        root_mult_array_f_x = [...
            0.7526  2;
            0.1236  3;
            1.1000  4;
            ];
        
        root_mult_array_f_y = [...
            0.1459  2;
            0.9876  1;
            ];
        
        m = 15;
    case '4'
        root_mult_array_f_x = ...
            [
            0.1564  5
            1.1256  2
            ];
        root_mult_array_f_y = ...
            [
            2.4456  1
            1.0576  3
            ];
        m = 11;
    otherwise
        error('Not a valid example number')
end

% Each row of polys x is a simple polynomial in p(x)

fprintf('Roots with respect to x:\n\n')
if length(root_mult_array_f_x) == 0 
    fprintf('None \n')
else
disp(root_mult_array_f_x)    
end

% Print the roots
fprintf('Roots with respect to y:\n\n')
if length(root_mult_array_f_y) == 0
    fprintf('None \n')
else
disp(root_mult_array_f_y)
end

% Print non-separable factors
fprintf('Non Separable factors \n\n')
if length(polys_xy) == 0
    fprintf('None \n')
else
    for i = 1:1:length(polys_xy)
    disp(polys_xy{i})
    end
end


polys_x = mult_roots_x(root_mult_array_f_x);
polys_y = mult_roots_y(root_mult_array_f_y);


polys_f = [polys_x; polys_y; polys_xy];

%PrintRoots_Bivariate(roots_f_x,roots_f_y,polys_xy);

fxy = BuildPoly_NonSeparable(polys_f);


end

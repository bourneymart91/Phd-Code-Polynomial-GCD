function [] = o_Intersection_Implicit_Parametric(ex_num,bool_preproc,low_rank_approx_method)
% Given an implicitly defined curve C1 := f(x,y) = 0, and a parametrically 
% defined curve C2, compute the points of intersection.



t = sym('t');
x = sym('x');
y = sym('y');

% 
SetGlobalVariables(bool_preproc,low_rank_approx_method);



switch ex_num
    case '1'
        % Intersection of y = x^2+1 (as a paremtric curve) and y = 32x -47
        % defined implicitly.
        f_xy = ...
            [
            -47 -1;
            32  0
            ];
        
        
        
        gx_t = ...
            [
            0;
            1
            ];
        
        gy_t = ...
            [
            1
            0
            0
            0
            0   
            0
            1
            ];
        gw_t = ...
            [
            1
            ];
    case '2' 
        % Intersection of two circles, one defined implicitly and one
        % defined as a rational parametric.
        
        % fxy = *(x-2)^2 + y^2 - 1 = 0
        f_xy = ...
            [
            4 0 -1;
            -4  0 0;
            1 0 0;
            ];
        
        gx_t = ...
            [
            1;
            0;
            -1;
            ];
        gy_t = ...
            [
            0;
            2;
            ];
        gw_t = ...
            [
            1;
            0;
            1;
            ];
        
end

% Get Symbolic expression for gx_t 
n1 = GetDegree(gx_t);
symbolic_gx_t = sum(diag(t.^(0:1:n1)) * gx_t);

% Get symbolic expression for gy_t
n2 = GetDegree(gy_t);
symbolic_gy_t = sum(diag(t.^(0:1:n2)) * gy_t);

% Get Symbolic expression for gw_t
n3 = GetDegree(gw_t);
symbolic_gw_t = sum(diag(t.^(0:1:n3)) * gw_t);

[m1,m2] = GetDegree_Bivariate(f_xy);

symbolic_fxy = diag(x.^(0:1:m1)) * f_xy * diag(y.^(0:1:m2));
symbolic_c3 = expand(sum(sum(subs(symbolic_fxy,{x,y}, ...
    {(symbolic_gx_t/symbolic_gw_t),(symbolic_gy_t./symbolic_gw_t)}))));

% Get as vector of coefficients of C3
c3 = fliplr(sym2poly(symbolic_c3))';

% Get the matrix of the roots and multiplicities
root_mult_mat = o_roots_mymethod(c3)

[nRoots,~] = size(root_mult_mat);

% For each root t 
for i = 1:1:nRoots
    
   % Get the root
   root = root_mult_mat(i);
   % Evaluate x(t) and y(t)
   x = polyval(flipud(gx_t),root);
   y = polyval(flipud(gy_t),root);
   coordinate = [x,y]
end



end
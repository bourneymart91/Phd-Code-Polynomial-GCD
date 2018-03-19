
function [fxy_matrix_exact,gxy_matrix_exact,...
    uxy_matrix_exact, vxy_matrix_exact,...
    dxy_matrix_exact,...
    m,m1,m2,...
    n,n1,n2,...
    t,t1,t2] = Examples_GCD_FromRoots(ex)

% NOTE -THESE EXAMPLES ASSUME SMALLEST POWER FIRST
f_roots_x = [];
f_roots_y = [];
f_roots_xy = [];
g_roots_x = [];
g_roots_y = [];
g_roots_xy = [];
d_roots_x = [];
d_roots_y = [];
d_roots_xy = [];
u_roots_x = [];
u_roots_y = [];
u_roots_xy = [];
v_roots_x = [];
v_roots_y = [];
v_roots_xy = [];

switch ex
    
    case 'template'
        f_roots_x = ...
            [
            ];
        f_roots_y = ...
            [
            ];
        
        f_roots_xy = ...
            [
            ];
        
        g_roots_x = ...
            [
            ];
        
        g_roots_y = ...
            [
            ];
        
        g_roots_xy = ...
            [
            ];
        d_roots_x = ...
            [
            ];
        d_roots_y = ...
            [
            ];
        d_roots_xy = ...
            [
            ];
        u_roots_x = ...
            [
            ];
        u_roots_y = ...
            [
            ];
        
        
        v_roots_x = ...
            [
            ];
        
        v_roots_y = ...
            [
            ];
        
    case '1'
        f_roots_x = ...
            [
            -0.5161 5
            -7.1052 5
            -0.1132 3
            ];
        f_roots_y = ...
            [
            0.7725  3
            ];
        
        f_roots_xy = ...
            [
            ];
        
        g_roots_x = ...
            [
            -0.5161 5;
            -7.1052 5;
            -2.0476 7;
            -8.8614 7;
            
            ];
        
        g_roots_y = ...
            [
            0.7725  3;
            ];
        
        g_roots_xy = ...
            [
            ];
        d_roots_x = ...
            [
            -0.5161 5
            -7.1052 5
            ];
        d_roots_y = ...
            [
            0.7725  3
            ];
        d_roots_xy = ...
            [
            ];
        u_roots_x = ...
            [
            -0.1132 3
            ];
        u_roots_y = ...
            [
            ];
        
        
        v_roots_x = ...
            [
            -2.0476     7;
            -8.8614     7;
            ];
        
        v_roots_y = ...
            [
            ];
        
        m = 16;
        n = 27;
        t = 13;
        
    case '2'
        f_roots_x = ...
            [
            1.456   1
            0.567   2
            0.927   3
            ];
        
        f_roots_y = ...
            [
            ];
        
        f_roots_xy = ...
            [
            ];
        
        g_roots_x = ...
            [
            1.456   1
            0.567   2
            0.427   3
            ];
        
        g_roots_y = ...
            [
            ];
        
        
        
        
        
        d_roots_x = ...
            [
            1.456   1
            0.567   2
            ];
        
        d_roots_y = ...
            [
            ];
        
        d_roots_xy = [];
        
        u_roots_x = ...
            [
            0.927   3
            ];
        
        u_roots_y = ...
            [
            ];
        
        u_roots_xy = [];
        
        v_roots_x = ...
            [
            0.427   3
            ];
        
        v_roots_y = ...
            [
            ];
        
        v_roots_xy = [];
        
        m = 6;
        n = 6;
        t = 3;
    case '3'
        f_roots_x = ...
            [
            1.456   1
            0.567   2
            0.927   3
            ];
        
        f_roots_y = ...
            [
            1.000   1
            1.192   2
            1.752   2
            ];
        
        f_roots_xy = ...
            [
            0.1     1;
            1       0
            ];
        
        g_roots_x = ...
            [
            1.456   1
            0.567   2
            0.427   3
            ];
        
        g_roots_y = ...
            [
            1.927   1
            1.192   2
            ];
        
        g_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{2,1} =...
            [
            1.7     1;
            1       0;
            ]
        
        
        
        d_roots_x = ...
            [
            1.456   1
            0.567   2
            ];
        
        d_roots_y = ...
            [
            1.192   2;
            ];
        
        d_roots_xy = ...
            [
            0.1     1;
            1       0
            ];
        
        u_roots_x = ...
            [
            0.927   3
            ];
        
        u_roots_y = ...
            [
            1.000   1;
            1.752   2;
            ];
        
        u_roots_xy = [];
        
        v_roots_x = ...
            [
            0.427   3
            ];
        
        v_roots_y = ...
            [
            1.927   1
            ];
        
        v_roots_xy = ...
            [
            1.7     1;
            1       0;
            ];
        
        m = 12;
        n = 11;
        t = 6;
        
    case '4'
        f_roots_x = ...
            [
            1.4   1
            0.5   2
            0.9   3
            ];
        
        f_roots_y = ...
            [
            1.3   1
            1.1   2
            1.7   2
            ];
        
        f_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        f_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        f_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_x = ...
            [
            1.4   1
            0.5   2
            0.4   3
            ];
        
        g_roots_y = ...
            [
            1.9   1
            1.1   2
            ];
        
        g_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{4,1} =...
            [
            1.7     1;
            1       0;
            ];
        
        
        
        d_roots_x = ...
            [
            1.4   1
            0.5   2
            ];
        
        d_roots_y = ...
            [
            1.1   2;
            ];
        
        d_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        
        u_roots_x = ...
            [
            0.927   3
            ];
        
        u_roots_y = ...
            [
            1.3   1;
            1.7   2;
            ];
        
        u_roots_xy = [];
        
        v_roots_x = ...
            [
            0.4   3
            ];
        
        v_roots_y = ...
            [
            1.9   1
            ];
        
        v_roots_xy = ...
            [
            1.7     1;
            1       0;
            ];
        
        m = 14;
        n = 13;
        t = 8;
        
    case '5'
        f_roots_x = ...
            [
            1.456   1
            0.567   2
            0.927   3
            ];
        
        f_roots_y = ...
            [
            1.000   1
            1.192   2
            1.752   2
            ];
        
        f_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        f_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        f_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        
        g_roots_x = ...
            [
            1.456   1
            0.567   2
            0.427   3
            ];
        
        g_roots_y = ...
            [
            1.927   1
            1.192   2
            ];
        
        g_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{4,1} =...
            [
            1.7     1;
            1       0;
            ];
        
        
        
        d_roots_x = ...
            [
            1.456   1
            0.567   2
            ];
        
        d_roots_y = ...
            [
            1.192   2;
            ];
        
        d_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        
        u_roots_x = ...
            [
            0.927   3
            ];
        
        u_roots_y = ...
            [
            1.000   1;
            1.752   2;
            ];
        
        u_roots_xy = [];
        
        v_roots_x = ...
            [
            0.427   3
            ];
        
        v_roots_y = ...
            [
            1.927   1
            ];
        
        v_roots_xy = ...
            [
            1.7     1;
            1       0;
            ];
        
        m = 14;
        n = 13;
        t = 8;
        
    case '6'
        f_roots_x = ...
            [
            1.456   1
            0.567   2
            0.927   3
            ];
        
        f_roots_y = ...
            [
            1.000   1
            1.192   2
            1.752   2
            ];
        
        f_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        f_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        f_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        
        g_roots_x = ...
            [
            1.456   1
            0.567   2
            0.427   2
            ];
        
        g_roots_y = ...
            [
            1.927   2
            1.192   2
            ];
        
        g_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{4,1} =...
            [
            1.7     1;
            1       0;
            ];
        
        
        
        d_roots_x = ...
            [
            1.456   1
            0.567   2
            ];
        
        d_roots_y = ...
            [
            1.192   2;
            ];
        
        d_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        
        u_roots_x = ...
            [
            0.927   3
            ];
        
        u_roots_y = ...
            [
            1.000   1;
            1.752   2;
            ];
        
        u_roots_xy = [];
        
        v_roots_x = ...
            [
            0.427   2
            ];
        
        v_roots_y = ...
            [
            1.927   2
            ];
        
        v_roots_xy = ...
            [
            1.7     1;
            1       0;
            ];
        
        m = 14;
        n = 13;
        t = 8;
        
    case '7'
        f_roots_x = ...
            [
            1.456   1;
            0.567   2;
            0.927   3;
            ];
        
        f_roots_y = ...
            [
            1.000   1;
            1.192   2;
            1.752   2;
            ];
        
        f_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0;
            ];
        f_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0;
            ];
        f_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0;
            ];
        
        g_roots_x = ...
            [
            1.456   1;
            0.567   2;
            0.427   3;
            ];
        
        g_roots_y = ...
            [
            1.927   4;
            1.192   2
            ];
        
        
        
        g_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{4,1} =...
            [
            1.7     1;
            1       0;
            ];
        
        
        
        d_roots_x = ...
            [
            1.456   1
            0.567   2
            ];
        
        d_roots_y = ...
            [
            1.192   2;
            ];
        
        d_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        
        u_roots_x = ...
            [
            0.927   3
            ];
        
        u_roots_y = ...
            [
            1.000   1;
            1.752   2;
            ];
        
        u_roots_xy = [];
        
        v_roots_x = ...
            [
            0.427   3
            ];
        
        v_roots_y = ...
            [
            1.927   4
            ];
        
        v_roots_xy = ...
            [
            1.7     1;
            1       0;
            ];
        
        m = 14;
        n = 13;
        t = 8;
        
    case '8'
        
        f_roots_x = ...
            [
            1.4   1
            0.5   2
            0.9   3
            ];
        
        f_roots_y = ...
            [
            1.0   1
            1.1   2
            1.7   2
            ];
        
        f_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        f_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        f_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        
        g_roots_x = ...
            [
            1.4   1
            0.5   2
            0.4   3
            ];
        
        g_roots_y = ...
            [
            1.9   1
            1.1   2
            ];
        
        g_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{4,1} =...
            [
            1.7     1;
            1       0;
            ];
        
        
        
        d_roots_x = ...
            [
            1.4   1
            0.5   2
            ];
        
        d_roots_y = ...
            [
            1.1   2;
            ];
        
        d_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        
        u_roots_x = ...
            [
            0.9   3
            ];
        
        u_roots_y = ...
            [
            1.0   1;
            1.7   2;
            ];
        
        u_roots_xy = [];
        
        v_roots_x = ...
            [
            0.4   3
            ];
        
        v_roots_y = ...
            [
            1.9   1
            ];
        
        v_roots_xy = ...
            [
            1.7     1;
            1       0;
            ];
        
        m = 14;
        n = 13;
        t = 8;
    case '9'
        f_roots_x = ...
            [
            1.456   1
            0.567   2
            0.927   3
            ];
        
        f_roots_y = ...
            [
            1.000   1
            1.192   2
            1.752   2
            ];
        
        f_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        f_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        f_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        
        
        g_roots_x = ...
            [
            1.456   1
            0.567   2
            0.427   3
            ];
        
        g_roots_y = ...
            [
            1.927   1
            1.192   2
            ];
        
        g_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{4,1} =...
            [
            1.7     1;
            1       0;
            ];
        
        
        
        d_roots_x = ...
            [
            1.456   1
            0.567   2
            ];
        
        d_roots_y = ...
            [
            1.192   2;
            ];
        
        d_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        
        u_roots_x = ...
            [
            0.927   3
            ];
        
        u_roots_y = ...
            [
            1.000   1;
            1.752   2;
            ];
        
        u_roots_xy = [];
        
        v_roots_x = ...
            [
            0.427   3
            ];
        
        v_roots_y = ...
            [
            1.927   1
            ];
        
        v_roots_xy = ...
            [
            1.7     1;
            1       0;
            ];
        
        m = 15;
        n = 14;
        t = 9;
    case '10'
        f_roots_x = ...
            [
            1.456   1
            0.567   2
            0.927   3
            ];
        
        f_roots_y = ...
            [
            1.000   1
            1.192   2
            1.752   2
            ];
        
        f_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        f_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        f_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        f_roots_xy{4,1} = ...
            [
            0.5     1;
            1       0;
            ];
        
        
        g_roots_x = ...
            [
            1.456   1
            0.567   2
            0.427   3
            ];
        
        g_roots_y = ...
            [
            1.927   1
            1.192   2
            ];
        
        g_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{4,1} =...
            [
            1.7     1;
            1       0;
            ];
        
        
        
        d_roots_x = ...
            [
            1.456   1
            0.567   2
            ];
        
        d_roots_y = ...
            [
            1.192   2;
            ];
        
        d_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        
        u_roots_x = ...
            [
            0.927   3
            ];
        
        u_roots_y = ...
            [
            1.000   1;
            1.752   2;
            ];
        
        u_roots_xy = ...
            [
            0.5     1;
            1       0;
            ];
        
        v_roots_x = ...
            [
            0.427   3
            ];
        
        v_roots_y = ...
            [
            1.927   1
            ];
        
        v_roots_xy = ...
            [
            1.7     1;
            1       0;
            ];
        
        m = 16;
        n = 14;
        t = 9;
        
    case '11'
        
        % Roots of f with respect to x
        f_roots_x = ...
            [
            2.4   1
            1.5   2
            1.9   3
            ];
        % Roots of f with respect to y
        f_roots_y = ...
            [
            0.3   1
            1.1   2
            1.7   2
            ];
        % Roots of f with respect to xy
        f_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0;
            ];
        
        f_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0;
            ];
        f_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0;
            ];
        % Roots of g with respect to x
        g_roots_x = ...
            [
            2.4   1
            1.5   2
            2.9   3
            ];
        
        % Roots of g with respect to y
        g_roots_y = ...
            [
            1.9   1
            1.1   2
            ];
        
        % Roots of f with respect to xy
        g_roots_xy{1,1} = ...
            [
            0.1     1;
            1     0
            ];
        g_roots_xy{2,1} = ...
            [
            0.1     1;
            1     0
            ];
        g_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        g_roots_xy{4,1} =...
            [
            1.7     1;
            1     0;
            ];
        
        
        % Roots of d with respect to x
        d_roots_x = ...
            [
            2.4   1
            1.5   2
            ];
        % Roots of d with respect to y
        d_roots_y = ...
            [
            1.1   2;
            ];
        
        % Roots of d with respect to xy
        d_roots_xy{1,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{2,1} = ...
            [
            0.1     1;
            1       0
            ];
        d_roots_xy{3,1} = ...
            [
            0.1     1;
            1       0
            ];
        
        
        % Roots of u with respect to x
        u_roots_x = ...
            [
            1.9   3
            ];
        
        % Roots of u with respect to y
        u_roots_y = ...
            [
            0.3   1;
            1.7   2;
            ];
        
        
        % Roots of v with respect to x
        v_roots_x = ...
            [
            2.9   3
            ];
        
        % Roots of v with respect to y
        v_roots_y = ...
            [
            1.9   1
            ];
        
        % Roots of v with respect to xy
        v_roots_xy = ...
            [
            1.7     1;
            1       0;
            ];
        
        m = 14;
        n = 13;
        t = 8;
        
    case '12'
        f_roots_x =  ...
            [
            0.1 2
            0.2 1
            0.5 1
            ];
        g_roots_x = ...
            [
            0.1 2
            0.3 1
            0.5 1
            ];
        d_roots_x = ...
            [
            0.1 2
            0.5 1
            ];
        u_roots_x = ...
            [
            0.2 1
            ];
        v_roots_x = ...
            [
            0.3 1
            ];
        m = 4;
        n = 4;
        t = 3;
        
   
    otherwise
        error('Not a valid example number')
        
end

f_roots_x = mult_roots_x(f_roots_x);
f_roots_y = mult_roots_y(f_roots_y);
g_roots_x = mult_roots_x(g_roots_x);
g_roots_y = mult_roots_y(g_roots_y);
d_roots_x = mult_roots_x(d_roots_x);
d_roots_y = mult_roots_y(d_roots_y);
u_roots_x = mult_roots_x(u_roots_x);
u_roots_y = mult_roots_y(u_roots_y);
v_roots_x = mult_roots_x(v_roots_x);
v_roots_y = mult_roots_y(v_roots_y);

f_roots = [f_roots_x; f_roots_y ; f_roots_xy];
g_roots = [g_roots_x; g_roots_y ; g_roots_xy];
d_roots = [d_roots_x; d_roots_y ; d_roots_xy];
u_roots = [u_roots_x; u_roots_y ; u_roots_xy];
v_roots = [v_roots_x; v_roots_y ; v_roots_xy];

fxy_matrix_exact = BuildPoly_NonSeparable(f_roots);
gxy_matrix_exact = BuildPoly_NonSeparable(g_roots);
dxy_matrix_exact = BuildPoly_NonSeparable(d_roots);
uxy_matrix_exact = BuildPoly_NonSeparable(u_roots);
vxy_matrix_exact = BuildPoly_NonSeparable(v_roots);

[r,c] = size(fxy_matrix_exact);
m1 = r -1;
m2 = c -1;

[r,c] = size(gxy_matrix_exact);
n1 = r -1;
n2 = c -1;

[r,c] = size(dxy_matrix_exact);
t1 = r-1;
t2 = c-1;


end

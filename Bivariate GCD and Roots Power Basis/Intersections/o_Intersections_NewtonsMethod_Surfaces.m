function [] = o_Intersections_NewtonsMethod_Surfaces()
% Compute the roots of a bivariate polynomial f(x,y) by Newtons Method


% Intial approximations of x,y
ex_num = '1';
x0(1) = 0.8;
y0(1) = 0.5;
z0(1) = 0.5;

syms x y z

switch ex_num
    case '1'
        f_sym = x^2 + y^2 + z^2 - 10 ;
    
        g_sym = x^2 + y^2 + (z - 5)^2 - 9 ;
end

test = figure();
hold on
ezimplot3(f_sym)
ezimplot3(g_sym)
hold off

ite = 1;



f_eval = subs(f_sym,[x,y,z],[x0(ite),y0(ite),z0(ite)]);
g_eval = subs(g_sym,[x,y,z],[x0(ite),y0(ite),z0(ite)]);

cond1(ite) = abs(f_eval) - abs(g_eval);
cond2(ite) = abs(f_eval);
cond3(ite) = abs(g_eval);


MAX_ERROR = 1E-14;
MAX_ITERATIONS = 20;

while ((cond1(ite) > MAX_ERROR) || (cond2(ite) > MAX_ERROR) || (cond3(ite) > MAX_ERROR)) && (ite < MAX_ITERATIONS)
    
    % evaluate f1 at original guess
    f_eval = subs(f_sym,[x,y,z],[x0(ite),y0(ite),z0(ite)]);
    g_eval = subs(g_sym,[x,y,z],[x0(ite),y0(ite),z0(ite)]);
    
    % Get partial derivatives
    df_dx_sym = diff(f_sym,x);
    df_dy_sym = diff(f_sym,y);
    df_dz_sym = diff(f_sym,z);
    
    dg_dx_sym = diff(g_sym,x);
    dg_dy_sym = diff(g_sym,y);
    dg_dz_sym = diff(g_sym,z);
    
    df_dx_eval = double(subs(df_dx_sym,[x,y,z],[x0(ite),y0(ite),z0(ite)]));
    df_dy_eval = double(subs(df_dy_sym,[x,y,z],[x0(ite),y0(ite),z0(ite)]));
    df_dz_eval = double(subs(df_dz_sym,[x,y,z],[x0(ite),y0(ite),z0(ite)]));
    
    dg_dx_eval = double(subs(dg_dx_sym,[x,y,z],[x0(ite),y0(ite),z0(ite)]));
    dg_dy_eval = double(subs(dg_dy_sym,[x,y,z],[x0(ite),y0(ite),z0(ite)]));
    dg_dz_eval = double(subs(dg_dz_sym,[x,y,z],[x0(ite),y0(ite),z0(ite)]));
    
    mat = [...
        df_dx_eval df_dy_eval df_dz_eval;
        dg_dx_eval dg_dy_eval dg_dz_eval
        ];
    
    vec = -1.* [f_eval ; g_eval];
    
    solution_vec = pinv(mat)*vec;
    
    % Get the small changes in x_{1} and x_{2}
    delta_x = solution_vec(1);
    delta_y = solution_vec(2);
    delta_z = solution_vec(3);
    
    % add the small changes to x_{1} and x_{2}
    x0(ite+1) = x0(ite) + delta_x;
    y0(ite+1) = y0(ite) + delta_y;
    z0(ite+1) = z0(ite) + delta_z;
    
    % Increment iteration counter
    ite = ite + 1;
    
    cond1(ite) = abs(f_eval) - abs(g_eval);
    cond2(ite) = abs(f_eval);
    cond3(ite) = abs(g_eval);
end

figure()
hold on
plot(log10(cond1),'-s','DisplayName','Condition 1');
plot(log10(cond2),'-s','DisplayName','Condition 2');
plot(log10(cond3),'-s','DisplayName','Condition 3');
xlabel('Iteration')
ylabel('log_{10}')
legend(gca,'show')
hold off

figure()
hold on
plot(log10(f_eval),'-s','DisplayName','Evaluation of f at (x_{i},y_{i})')
plot(log10(g_eval),'-s','DisplayName','Evaluation of g at (x_{i},y_{i})')
xlabel('Iteration number')
ylabel('log_{10}f(x,y) - g(x,y)')
legend(gca,'show')
hold off

figure(test)
hold on
plot3(x0,y0,z0,'-s','DisplayName','Approximation of intersection point')
xlabel('x');
ylabel('y');
zlabel('z');
title('Approximate values of (x,y) by Newton root finding')
hold off

figure()
hold on
plot3(log10(x0),log10(y0),log10(z0),'-s')
xlabel('log_{10}(x)');
ylabel('log_{10}(y)');
zlabel('log_{10}(z)');
title('Approximate values of log(x,y) by Newton root finding')
hold off

format long
display(cond1)
display(cond2)
display(cond3)
display(x0(ite))
display(y0(ite))

end



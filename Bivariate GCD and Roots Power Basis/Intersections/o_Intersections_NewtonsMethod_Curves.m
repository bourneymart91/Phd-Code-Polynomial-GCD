function [] = o_Intersections_NewtonsMethod_Curves()
% Given an initial approximation of an intersection point, Compute the
% intersection point to within a specified threshold by Newtons Method,
% Extended from Newton's Root finding method


syms x y

ex_num = '3';

switch ex_num
    case '1'
        f_sym = (x)^2 + y^2 - 10;
        g_sym = (x-1)^2 + y^2 -8;
        
    case '2'
        f_sym = x^2 + y^2 - 4;
        g_sym = x^2 + y^2 - 6*x + 4;  
                      
    case '3'
        f_sym = x^2 + y^2 -4 ;
        g_sym = y -0 +0*x;
    otherwise 
        error('error')
end


test = figure();
hold on
ezplot(f_sym)
ezplot(g_sym)
hold off

ite = 1;

x0(ite) = 1;
y0(ite) = -1;

f_eval = subs(f_sym,[x,y],[x0(ite),y0(ite)]);
g_eval = subs(g_sym,[x,y],[x0(ite),y0(ite)]);

cond1(ite) = abs(f_eval) - abs(g_eval);
cond2(ite) = abs(f_eval);
cond3(ite) = abs(g_eval);

MAX_ERROR = 1E-14;
MAX_ITERATIONS = 20;

while ((cond1(ite) > MAX_ERROR) || (cond2(ite) > MAX_ERROR) || (cond3(ite) > MAX_ERROR)) && (ite < MAX_ITERATIONS)
    
    % evaluate f1 at original guess
    f_eval = double(subs(f_sym,[x,y],[x0(ite),y0(ite)]));
    g_eval = double(subs(g_sym,[x,y],[x0(ite),y0(ite)]));
    
    % Get partial derivatives
    df_dx_sym = diff(f_sym,x);
    df_dy_sym = diff(f_sym,y);
    dg_dx_sym = diff(g_sym,x);
    dg_dy_sym = diff(g_sym,y);
        
    df_dx_eval = double(subs(df_dx_sym,[x,y],[x0(ite),y0(ite)]));
    df_dy_eval = double(subs(df_dy_sym,[x,y],[x0(ite),y0(ite)]));
    dg_dx_eval = double(subs(dg_dx_sym,[x,y],[x0(ite),y0(ite)]));
    dg_dy_eval = double(subs(dg_dy_sym,[x,y],[x0(ite),y0(ite)]));
    
    mat = [...
        df_dx_eval df_dy_eval;
        dg_dx_eval dg_dy_eval
        ];
    
    vec = -1.* [f_eval ; g_eval];
    
    solution_vec = pinv(mat)*vec;
    
    % Get the small changes in x_{1} and x_{2}
    delta_x1 = solution_vec(1);
    delta_x2 = solution_vec(2);
    
    % add the small changes to x_{1} and x_{2}
    x0(ite+1) = x0(ite) + delta_x1;
    y0(ite+1) = y0(ite) + delta_x2;
    
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
plot(x0,y0,'-s','DisplayName','Approximation of intersection point')
xlabel('x');
ylabel('y');
title('Approximate values of (x,y) by Newton root finding')
hold off

figure()
hold on
plot(log10(x0),log10(y0),'-s')
xlabel('log_{10}(x)');
ylabel('log_{10}(y)');
title('Approximate values of log(x,y) by Newton root finding')
hold off

format long
display(cond1)
display(cond2)
display(cond3)
display(x0(ite))
display(y0(ite))

end



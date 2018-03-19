function [] = o_Roots_NewtonsMethod()
% Compute the roots of a bivariate polynomial f(x,y) by Newtons Method


% Intial approximations of x,y
ex_num = '1';
x0(1) = 0.8;
y0(1) = 0.5;

syms x y

switch ex_num
    case '1'
        f_sym = x^2 + y^2 ;
    case '2'
        f_sym = x^3 + 3*x^2 + 2*x + 1;
end

df_dx_sym = diff(f_sym,x);
df_dy_sym = diff(f_sym,y);

% evaluate f(x,y) at x0 y0
ite = 1;
f_evaluated(ite) = subs(f_sym,[x,y],[x0,y0]);

% Get partial derivatives

df_dx_evaluated = subs(df_dx_sym,[x,y],[x0,y0]);
df_dy_evaluated = subs(df_dy_sym,[x,y],[x0,y0]);

% build jacobi
mat = [df_dx_evaluated df_dy_evaluated];

delta_x = 0;
delta_y = 0;



MAX_ITERATIONS = 50;
MAX_ERROR = 1e-12;

while abs(f_evaluated(ite)) > MAX_ERROR && ite < MAX_ITERATIONS

    x0(ite+1) = x0(ite) + delta_x;
    y0(ite+1) = y0(ite) + delta_y;

    df_dx_evaluated = subs(df_dx_sym,[x,y],[x0(ite+1),y0(ite+1)]);
    df_dy_evaluated = subs(df_dy_sym,[x,y],[x0(ite+1),y0(ite+1)]);

    f_evaluated(ite+1) = subs(f_sym,[x,y],[x0(ite+1),y0(ite+1)]);

    mat = [df_dx_evaluated df_dy_evaluated];


    delta = (pinv(mat) * (-1.*f_evaluated(ite+1)) );
    delta_x = delta(1);
    delta_y = delta(2);

    ite = ite + 1;

end

fprintf('End \n')
fprintf('Numer of iterations Required : %i \n',ite)
display(x0(ite))
display(y0(ite))

% Plot the evaluated f(x,y) for each (x_{i},y_{i})
figure_name = sprintf([mfilename ' : ' 'Evaluated f(x,y)']);
figure('Name',figure_name)
hold on
plot(log10(f_evaluated),'-s','DisplayName','f evaluated at (x_{i},y_{i})')
xlabel('Iteration number i')
ylabel('log_{10}f(x_{i},y_{i})')
hold off

% Plot (x_{i},y_{i}) for every iteration
figure_name = sprintf([mfilename ' : ' 'Values x_{i},y_{i}']);
figure('Name',figure_name)
hold on
plot(x0,y0,'-s')
xlabel('x');
ylabel('y');
title('Approximate values of (x,y) by Newton root finding')
hold off

figure()
hold on
scatter(log10(x0),log10(y0))
xlabel('log_{10}(x)');
ylabel('log_{10}(y)');
title('Approximate values of (x,y) by Newton root finding')
hold off


end
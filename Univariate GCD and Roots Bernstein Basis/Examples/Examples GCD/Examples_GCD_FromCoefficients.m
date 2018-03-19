function [fx,gx,dx,ux,vx] = Examples_GCD_FromCoefficients(ex_num)
%
% % Inputs
%
% ex_num : Example number
%
% % Outputs
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x)
%
% dx : (Vector) Coefficients of the GCD d(x)
%
% ux : (Vector) Coefficients of the polynomial u(x) given by f(x)/d(x)
%
% vx : (Vector) Coefficients of the polynomial v(x) given by g(x)/d(x)

addpath(genpath('../Examples'));

% Get roots and multiplicities from example file
[f_root_mult_arr, g_root_mult_arr, d_root_mult_arr, u_root_mult_arr, v_root_mult_arr] = ...
    GCD_Examples_Univariate_2Polys(ex_num);



PlotRoots(f_root_mult_arr, g_root_mult_arr, d_root_mult_arr);

% Get the coefficients of the polynomials f(x), g(x) and d(x).
fx = BuildPolyFromRootsSymbolic(f_root_mult_arr);
gx = BuildPolyFromRootsSymbolic(g_root_mult_arr);
dx = BuildPolyFromRootsSymbolic(d_root_mult_arr);
ux = BuildPolyFromRootsSymbolic(u_root_mult_arr);
vx = BuildPolyFromRootsSymbolic(v_root_mult_arr);

% Get the symbolic polynomials
fx_sym = GetSymbolicPolyFromSymbolicRoots(f_root_mult_arr);
gx_sym = GetSymbolicPolyFromSymbolicRoots(g_root_mult_arr);
dx_sym = GetSymbolicPolyFromSymbolicRoots(d_root_mult_arr);
ux_sym = GetSymbolicPolyFromSymbolicRoots(u_root_mult_arr);
vx_sym = GetSymbolicPolyFromSymbolicRoots(v_root_mult_arr);

display(fx_sym)
display(gx_sym)
display(dx_sym)
display(ux_sym)
display(vx_sym)

m = feval(symengine, 'degree', fx_sym);
n = feval(symengine, 'degree', gx_sym);
t = feval(symengine, 'degree', dx_sym);

fprintf('Degree f(x) : %i \n', m)
fprintf('Degree g(x) : %i \n', n)
fprintf('Degree d(x) : %i \n', t)

end

function PlotRoots(f_root_mult_arr, g_root_mult_arr, d_root_mult_arr)


figure()
hold on

const = 0.1;
% for each root
for i = 1 : 1: size(f_root_mult_arr, 1)
   
    color = 'r';
    factor = f_root_mult_arr(i,1);
    vCoeffs = coeffs(factor);
    root = vCoeffs(1);
    mult = f_root_mult_arr(i,2);
    circle(root, 0 ,mult * const, color);
    
    
end

for i = 1 : 1: size(g_root_mult_arr, 1)
   
    color = 'b';
    factor = g_root_mult_arr(i,1);
    vCoeffs = coeffs(factor);
    root = vCoeffs(1);
    mult = g_root_mult_arr(i,2);
    h = circle(root, 0 ,mult * const, color);
    
    
end

for i = 1 : 1: size(d_root_mult_arr, 1)
   
    color = 'g';
    factor = d_root_mult_arr(i,1);
    vCoeffs = coeffs(factor);
    root = vCoeffs(1);
    mult = d_root_mult_arr(i,2);
    h = circle(root, 0 ,mult * const, color);
    
    
end


grid on
hold off

end



function h = circle(x,y,r, color)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit, color);
    h2 = plot(x,y,'-s','color',color);
    hold off
end

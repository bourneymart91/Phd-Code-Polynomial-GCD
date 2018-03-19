function [] = Intersection_test()


close all; clc;

syms x

f_root_mult_arr = [...
    (x - 0.2) 4;
    (x + 0.1) 2
    ] ;
g_root_mult_arr = [...
    (x - 0.2)   4
    (x + 0.5)   1
    ];


fx = BuildPolyFromRootsSymbolic(f_root_mult_arr);
gx = BuildPolyFromRootsSymbolic(g_root_mult_arr);

CP_fx = GetControlPointsBernsteinPoly(fx);
CP_gx = GetControlPointsBernsteinPoly(gx);


% Plot from Bernstein form
Plot2DBezier({CP_fx, CP_gx});

% Plot from symbolic power basis form
fx_sym = GetSymbolicPolyFromSymbolicRoots(f_root_mult_arr);
gx_sym = GetSymbolicPolyFromSymbolicRoots(g_root_mult_arr);
PlotSymbolicPoly({fx_sym, gx_sym});




% Implicitize g(x)
CP_gx = GetControlPointsBernsteinPoly(gx);
[imp_gx, symbolic_gx] = ImplicitizeBezierBySylvester(CP_gx);
PlotSymbolicPoly({symbolic_gx});


% Substitute
hx = Substitute(CP_gx, fx);

% Get control points of h(x)
CP_hx = GetControlPointsBernsteinPoly(hx);

% Plot control points
Plot2DBezier({CP_hx});

matlab_roots = o_roots_Matlab(hx);
disp(matlab_roots);

o_roots_mymethod(hx);





sameaxes()

end




function PlotSymbolicPoly(arrPolys)

figure()
hold on
x_low = 0;
x_high = 1;

nPolys = size(arrPolys,2);

for i = 1 : 1 : nPolys
    fx = arrPolys{i};
    ezplot(fx, [x_low, x_high])
end

hold off

end



function CP = GetControlPointsBernsteinPoly(fx)

m = GetDegree(fx);
xVec = 0:(1/m):1;
CP = [xVec; fx'];

end
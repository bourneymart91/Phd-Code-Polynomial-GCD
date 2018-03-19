function [hx] = Deconvolve(fx, gx)
% DECONVOLVE Given two polynomails f(x,y) and g(x,y), computes the
% polynomail division f(x)/g(x) = h(x)
%
% % Inputs.
%
% fx : Coefficients of polynomial f(x)
%
% gx : Coefficients of polynomial g(x)
%
% % Outputs.
%
% hx : Coefficients of polynomial h(x)


% Get degree of polynomial f(x).
m = GetDegree(fx);

% Get degree of polynomial g(x).
n = GetDegree(gx);

% Build the matrix T_{m-n}(g(x))
C_g = BuildT1(gx, m-n);

% Solve T_{m-n}(g(x))*h = f
hx = SolveAx_b(C_g, fx);

end

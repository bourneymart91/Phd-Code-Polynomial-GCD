function h1 = Deconvolve(fx, gx)
% Given two polynomials, f(x) and g(x) of degrees m and n, perform
% deconvolution to obtain h(x) = f(x)/g(x), 
%
% Inputs.
%
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% gx : (Vector) Coefficients of polynomial g(x)
%
% Outputs.
%
%
% hx : (Vector) Coefficients of polynomial h(x)


% h1 : output polynomial = f(x)/g(x).

% Obtain degree of polynomial f(x)
m = GetDegree(fx);

% Get the degree of polynomial g(x)
n = GetDegree(gx);

% Build the matrix DCQ
DCQ = BuildDT1Q1(gx, m-n);

% Solve C(g) * h = f, for the unknown vector h.
h1 = SolveAx_b(DCQ, fx);

end
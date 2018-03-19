function hxy = Deconvolve_Bivariate(fxy, gxy)
% DECONVOLVE_BIVARIATE :
% Deconvolve the polynomials f(x,y) and g(x,y) to obtain h(x,y)
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
% 
% % Outputs
%
% hxy : (Matrix) Coefficients of polynomial h(x,y) = f(x,y) / g(x,y).


% Get degree of f(x,y) (Note: m1 = m2)
[m1, ~] = GetDegree_Bivariate(fxy);

% Get degree of g(x,y) (Note: n1 = n2)
[n1, ~] = GetDegree_Bivariate(gxy);

% Get m, the total degree of f(x,y)
m = m1;

% Get n, the total degree of g(x,y)
n = n1;

% %
% %
% Build the system Ax = b, and solve for x. Where A is the coefficient 
% matrix = D_{m}*T(g)*Q_{m-n}, and b is the vector of coefficients of
% f(x,y).

% Build the matrix D^{-1}_{n}
D = BuildD_2Polys(n, m - n);

% Build the matrix T_{m-n}(g)
T1 = BuildT1(gxy, n, m - n);

% Build the matrix Q_{m-n}
Q = BuildQ1(m - n);

% Build the matrix D^{-1}_{m}*T_{m-n}(g)Q_{m-n}
DTQ = D*T1*Q;

% %
% Build the rhs coefficient vector. Get vector of coefficients of f(x,y) 
% (Note - This includes trailing zeros from the lower right part of the 
% matrix, which must be deleted.

f = GetAsVector(fxy);

% Get number of coefficients in f(x,y)
nCoefficients_fxy = nchoosek(m+2, 2);

% Remove trailing zeros
v_fxy = f(1 : nCoefficients_fxy);

% % Solve for x
% Get coefficients of h(x,y) 
x_ls = SolveAx_b(DTQ, v_fxy);

% %
% (Note that x_ls only gives the upper left triangle of coefficients, so 
% append zeros to fill lower right triangle to form matrix of h(x,y))

% Get number of zeros in the coefficient matrix for polynomial h(x,y)
if (m - n + 1 > 1)
    nZeros_hxy = nchoosek(m - n + 1, 2);
else
    nZeros_hxy = 0;
end

% Append zeros to vector
v_hxy = ...
    [
    x_ls;
    zeros(nZeros_hxy, 1);
    ];

% Get the matrix of coefficients of h(x,y).
hxy = GetAsMatrix(v_hxy, m - n, m - n);



end


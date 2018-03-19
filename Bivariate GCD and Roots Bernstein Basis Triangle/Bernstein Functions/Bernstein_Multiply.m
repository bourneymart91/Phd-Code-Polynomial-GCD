function [hxy] = Bernstein_Multiply(fxy, gxy, m, n)
% Multiply the bivariate bernstein polynomials f(x,y) and g(x,y) of degrees
% m and n respectively.
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% m : (Int) Total degree of polynomial f(x,y)
%
% n : (Int) Total degree of polynomial g(x,y)

% Build the matrix D^{-1}_{m+n-k}
D = BuildD_2Polys(m,n);

% T_{n}(f)
T1 = BuildT1(fxy, m, n);

% Build the matrix Q_{n}
Q1 = BuildQ1(n);

% Get the matrix containing coefficients of g(x,y) as a vector.
g = GetAsVector(gxy);

% Get number of coefficients in g(x,y) = \nchoosek(n+2,2).
nCoeffs_gxy = nchoosek(n+2,2);

% Remove the zeros from the vector g
g = g(1:nCoeffs_gxy);

% Compute the vector h containing coefficients of h(x,y)
h = D*T1*Q1 * g;

% Append zeros to vector h so that the entries of h fill a matrix of size
% (m+n+1) x (m+n+1)
nZeros_hxy = nchoosek(m+n+1,2);
h = [h; zeros(nZeros_hxy,1)];

% Get the matrix of coefficients of h(x,y)
hxy = GetAsMatrix(h,m+n,m+n);

end
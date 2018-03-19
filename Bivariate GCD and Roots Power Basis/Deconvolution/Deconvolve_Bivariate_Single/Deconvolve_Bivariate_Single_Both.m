function [hxy_matrix] = Deconvolve_Bivariate_Single_Both(fxy, gxy, m, n)
% Perform polynomial deconvolution of the bivariate polynomials f(x,y) and 
% g(x,y) where the total degrees m and n of f and g respectively are known.
% This method is referred to as 'both' since BOTH the total degree and
% relative degrees of the polynomials are utilised.
%
% Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
%
% Outputs.
%
% hxy : (Matrix) Coefficients of polynomial h(x,y)

% Return the matrix of coefficients of the polynomial h, where h = f/g


% Get the degrees of polynomial f(x,y) and g(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

% Get the degrees of polynomial h(x,y)
degt_hxy = m - n;
deg1_hxy = m1 - n1;
deg2_hxy = m2 - n2;

% Get number of zeros in h(x,y)
nNonZeros_hxy = GetNumNonZeros(deg1_hxy,deg2_hxy,degt_hxy);
nZeros_hxy = (m1-n1 + 1) * (m2-n2 + 1) - nNonZeros_hxy;

% Get number of nonzeros in g*h = f(x,y)
nNonZeros_gh = GetNumNonZeros(m1, m2, m);

% % 
% %
% Build the matrix C(g)
C1 = BuildT1_Both_Bivariate(gxy, n, m-n, m1-n1, m2-n2);


% %
% %
% Create Right hand side vector f(x,y)

% Get the polynomial f(x,y) in vector form
f = GetAsVector_Version1(fxy);

% Remove zeros from f(x,y).
f = f(1:nNonZeros_gh, 1);

% Solve the Ax=b problem.
x_ls = SolveAx_b(C1, f);

% Get set of coefficients of h(x,y), including zeros to form matrix of
% coefficients
hxy_vec = ...
    [
    x_ls;
    zeros(nZeros_hxy,1);
    ];

% Get h(x,y) as a vector
hxy_matrix = GetAsMatrix_Version1(hxy_vec, m1-n1, m2-n2);


end


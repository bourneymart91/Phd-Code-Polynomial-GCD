function dxy = GetGCDCoefficients_Both_Bivariate_2Polys(fxy, gxy, uxy, vxy, m, n, k)
% Given the matrices of coefficients of f(x,y) and g(x,y), the quotient
% polynomials u(x,y) and v(x,y), and optimal values for alpha, theta_{1}
% and theta_{2}, calculate the coefficients of the GCD d(x,y).
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% uxy : (Matrix) Coefficients of polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of polynomial v(x,y)
%
% m n : (Int) (Int) Total degree of f(x,y) and g(x,y)
% 
% k : (Int) Total degree of d(x,y)
%
% % Outputs
%
% dxy : (Matrix) Coefficients of polynomial d(x,y)
%
% % Note : When both total and relative degree are used, always structure
% matrices and vectors in the 'Version 1' style.



% Calculate the GCD of two bivariate polynomials f(x,y) and g(x,y)

% %
% %
% Get Degree structures

% Get degrees of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

% Get degrees of polynomial u(x,y)
[m1_k1, m2_k2] = GetDegree_Bivariate(uxy);

% Get degrees of polynomial d(x,y)
t1 = m1 - (m1_k1);
t2 = m2 - (m2_k2);


% Build matrices C(u) and C(v)
% Get number of zeros in d(x,y)
nNonZeros_dxy = GetNumNonZeros(t1, t2, k);
nZeros_dxy = (t1+1) * (t2+1) - nNonZeros_dxy;

nNonZeros_fxy = GetNumNonZeros(m1, m2, m);
nNonZeros_gxy = GetNumNonZeros(n1, n2, n);

% Build the coefficient matrix
C1 = BuildT1_Both_Bivariate(uxy, m-k, k, t1, t2);
C2 = BuildT1_Both_Bivariate(vxy, n-k, k, t1, t2);

C = [C1 ; C2];

% % 
% % 
% Preprocess f(x,y) and g(x,y) and get in vector form

% Get fww_matrix as a vector
fxy_vec = GetAsVector_Version1(fxy);

% Remove the zeros from f(x,y) 
fxy_vec = fxy_vec(1:nNonZeros_fxy,:);

% get gww_matrix as a vector
gxy_vec = GetAsVector_Version1(gxy);

% Remove the zeros from g(x,y)
gxy_vec = gxy_vec(1:nNonZeros_gxy,:);

% Build the RHS vector
rhs_vec = [fxy_vec;
           gxy_vec];

% % Calculate the solution vector

% Calculate the x vector by pinv       
x = SolveAx_b(C,rhs_vec);


% Append the removed zeros
dxy_vec = ...
    [
        x;
        zeros(nZeros_dxy,1);
    ];

% Arrange d(x,y) into a matrix form based on its dimensions.
dxy = GetAsMatrix_Version1(dxy_vec, t1, t2);



end
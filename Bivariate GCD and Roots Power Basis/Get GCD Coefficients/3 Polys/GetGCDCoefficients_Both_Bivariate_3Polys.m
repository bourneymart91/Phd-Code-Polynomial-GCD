function dxy = GetGCDCoefficients_Both_Bivariate_3Polys(fxy, gxy, hxy, uxy, vxy, wxy, m, n, o, k)
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
% hxy : (Matrix) Coefficients of polynomial h(x,y)
% 
% uxy : (Matrix) Coefficients of polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of polynomial v(x,y)
%
% wxy : (Matrix) Coefficients of polynomial w(x,y)
%
% m : Total degree of f(x,y)
% 
% n : Total degree of g(x,y)
% 
% o : Total degree of w(x,y)
% 
% k : Total degree of d(x,y)
%
% % Outputs
%
% dxy : (Matrix) Coefficients of polynomial d(x,y)


% Calculate the GCD of two bivariate polynomials f(x,y) and g(x,y)

% %
% %
% Get Degree structures

% Get degrees of polynomial f(x,y), g(x,y) and h(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);

% Get degrees of polynomial u(x,y)
[m1_k1, m2_k2] = GetDegree_Bivariate(uxy);

% Get degrees of polynomial d(x,y)
t1 = m1-(m1_k1);
t2 = m2-(m2_k2);

% %
% %
% Build matrices C(u) and C(v)
% Get number of zeros in d(x,y)
nNonZeros_dxy = GetNumNonZeros(t1,t2,k);
nZeros_dxy = (t1+1) * (t2+1) - nNonZeros_dxy;

nNonZeros_ud = GetNumNonZeros(m1, m2, m);
nNonZeros_vd = GetNumNonZeros(n1, n2, n);
nNonZeros_wd = GetNumNonZeros(o1, o2, o);

C1 = BuildT1_Both_Bivariate(uxy, m-k, k, t1, t2);
C2 = BuildT1_Both_Bivariate(vxy, n-k, k, t1, t2);
C3 = BuildT1_Both_Bivariate(wxy, o-k, k, t1, t2);

C = [C1; C2; C3];

% % 
% % 
% Preprocess f(x,y) and g(x,y) and get in vector form

% Get fww_matrix as a vector and remove trailing zeros
fxy_vec = GetAsVector_Version1(fxy);
fxy_vec = fxy_vec(1:nNonZeros_ud,:);

% Get gww_matrix as a vector and remove trailing zeros
gxy_vec = GetAsVector_Version1(gxy);
gxy_vec = gxy_vec(1:nNonZeros_vd,:);

% Get hww_matrix as a vector and remove trailing zeros
hxy_vec = GetAsVector_Version1(hxy);
hxy_vec = hxy_vec(1:nNonZeros_wd,:);

% Build the RHS vector
rhs_vec = [fxy_vec;
           gxy_vec;
           hxy_vec];

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
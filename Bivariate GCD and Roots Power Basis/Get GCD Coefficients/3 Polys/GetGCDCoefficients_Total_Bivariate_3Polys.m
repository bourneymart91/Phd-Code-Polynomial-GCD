function dxy = GetGCDCoefficients_Total_Bivariate_3Polys(fxy, gxy, hxy, uxy, vxy, wxy, m, n, o, t)
% Given the matrices of coefficients of f(x,y) and g(x,y), the quotient
% polynomials u(x,y) and v(x,y), and optimal values for alpha, theta_{1}
% and theta_{2}, calculate the coefficients of the GCD d(x,y).
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of the polynomial h(x,y)
%
% uxy : (Matrix) Coefficients of the cofactor polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of the cofactor polynomial v(x,y)
%
% wxy : (Matrix) Coefficients of the cofactor polynomial w(x,y)
%
% m n o : (Int) (Int) (Int) : Total degree of f(x,y), g(x,y) and h(x,y)
%
% t : (Int) Total degree of polynomial d(x,y)

% % Build Matrix C
% Build the Cauchy matrix of coefficients of u(x,y), v(x,y) and w(x,y)
C1 = BuildT1_Total_Bivariate(uxy, m-t, t);
C2 = BuildT1_Total_Bivariate(vxy, n-t, t);
C3 = BuildT1_Total_Bivariate(wxy, o-t, t);

% Build the RHS vector of coefficients of f and g
C = [C1 ; C2 ;  C3];

% Get number of coefficients in f(x,y), g(x,y) and h(x,y) when expressed 
% in terms of total degree.
nCoefficients_fxy = nchoosek(m+2,2);
nCoefficients_gxy = nchoosek(n+2,2);
nCoefficients_hxy = nchoosek(o+2,2);


% % Build vector f(x,y)
% Pad f so that it is in terms of its total degree.
[m1, m2] = GetDegree_Bivariate(fxy);
padd = zeros(m+1,m+1);
padd(1:m1+1,1:m2+1) = fxy;
fxy = padd;

% % Build vector of coefficients of g(w,w)
% Pad g so that it is in terms of its total degree.
[n1, n2] = GetDegree_Bivariate(gxy);
padd = zeros(n+1,n+1);
padd(1:n1+1,1:n2+1) = gxy;
gxy = padd;

% % Build vector of coefficients of g(w,w)
% Pad g so that it is in terms of its total degree.
[o1, o2] = GetDegree_Bivariate(hxy);
padd = zeros(o+1,o+1);
padd(1:o1+1,1:o2+1) = hxy;
hxy = padd;

% Get fxy_matrix as a vector and remove the zeros associated with the 
% polynomial by total degree.
fxy_vec = GetAsVector_Version1(fxy);
fxy_vec = fxy_vec(1:nCoefficients_fxy);

% Get gww_matrix as a vector and remove trailing zeros
gxy_vec = GetAsVector_Version1(gxy);
gxy_vec = gxy_vec(1:nCoefficients_gxy);

% Get gxy_matrix as a vector and remove trailing zeros
hxy_vec = GetAsVector_Version1(hxy);
hxy_vec = hxy_vec(1:nCoefficients_hxy);

% % Build the RHS vector
rhs_vec = [...
    fxy_vec;...
    gxy_vec;...
    hxy_vec
    ];

% Calculate the vector x
x = SolveAx_b(C,rhs_vec);
dxy_vec = x;

% Get the residual associated with the solution x. (Small residual implies good approximation)
residual = pinv(C)*rhs_vec - x;

% Padd d(w,w) with zeros so that it can be put back into matrix form
dxy_vec = [dxy_vec ; zeros(nchoosek(t+2-1,2),1)];

% Arrange dw into a matrix form based on its dimensions.
dxy = GetAsMatrix_Version1(dxy_vec,t,t);


end
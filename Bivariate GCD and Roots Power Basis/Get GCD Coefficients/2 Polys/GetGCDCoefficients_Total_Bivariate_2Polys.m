function dxy = GetGCDCoefficients_Total_Bivariate_2Polys(fxy, gxy, uxy, vxy, m, n, t)
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
% uxy : (Matrix) Coefficients of the polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of the polynomial h(x,y)
%
% m : Total degree of polynomial f(x,y)
%
% n : Total degree of polynomial g(x,y)
%
% t : Total degree of polynomial d(x,y)
%
% % Outputs
%
% dxy : (Matrix) Coefficients of the polynomial d(x,y)

% % Build Matrix C
% Build the Cauchy matrix of coefficients of u(w,w)
C1 = BuildT1_Total_Bivariate(uxy, m-t, t);

% Build the Cauchy matrix of coefficients of v(w,w)
C2 = BuildT1_Total_Bivariate(vxy, n-t, t);

% Build the RHS vector of coefficients of f and g
C = [C1 ; C2];

% Get number of coefficients in f(x,y) in terms of total degree.
nCoefficients_fxy = nchoosek(m+2,2);

% Get number of coefficietns in g(x,y) in terms of total degree
nCoefficients_gxy = nchoosek(n+2,2);

% % Build vector f(x,y)
% Pad f so that it is in terms of its total degree.
[m1, m2] = GetDegree_Bivariate(fxy);
padd = zeros(m+1,m+1);
padd(1:m1+1, 1:m2+1) = fxy;
fxy = padd;

% % Build vector of coefficients of g(w,w)
% Pad g so that it is in terms of its total degree.
[n1, n2] = GetDegree_Bivariate(gxy);
padd = zeros(n+1,n+1);
padd(1:n1+1,1:n2+1) = gxy;
gxy = padd;

% Get fxy_matrix as a vector
fxy_vec = GetAsVector_Version1(fxy);

% Remove the zeros associated with the polynomial by total degree
fxy_vec = fxy_vec(1:nCoefficients_fxy);

% get gww_matrix as a vector
gxy_vec = GetAsVector_Version1(gxy);

% Remove the zeros associated with the polynomial by total degree
gxy_vec = gxy_vec(1:nCoefficients_gxy);

% % Build the RHS vector
rhs_vec = [...
    fxy_vec;...
    gxy_vec...
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
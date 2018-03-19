function dxy = GetGCDCoefficients_Relative_Bivariate_3Polys(fxy, gxy, hxy, uxy, vxy, wxy)
% Given the matrices of coefficients of f(x,y) and g(x,y), the quotient
% polynomials u(x,y) and v(x,y), and optimal values for alpha, theta_{1}
% and theta_{2}, calculate the coefficients of the GCD d(x,y), by forming
% the matrix-vector product C(u)*d = f and C(v)*d = g, and solve for d. 
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of the polynomial h(x,y)
%
% uxy : (Matrix) Coefficients of the polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of the polynomial v(x,y)
% 
% wxy : (Matrix) Coefficients of the polynomial w(x,y)
%
% % Outputs
%
% dxy :(Matrix) Coefficients of the polynomial d(x,y)



% Calculate the GCD of two bivariate polynomials f(x,y) and g(x,y)

% Get degree of polynomial f(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degree of polynomial u(x,y) with respect to x and y
[m1_t1, m2_t2] = GetDegree_Bivariate(uxy);

% Get degrees of polynomial d(x,y)
t1 = m1-(m1_t1);
t2 = m2-(m2_t2);

% %
% %
% Build Matrix C

% Build the Cauchy matrix of coefficients of u(\omega_{1},\omega_{2})
% Build the Cauchy matrix of coefficients of v(\omega_{1},\omega_{2})
% Build the Cauchy matrix of coefficients of w(\omega_{1},\omega_{2})
C1 = BuildT1_Relative_Bivariate_Version1(uxy, t1, t2);
C2 = BuildT1_Relative_Bivariate_Version1(vxy, t1, t2);
C3 = BuildT1_Relative_Bivariate_Version1(wxy, t1, t2);

% Build the RHS vector of coefficients of f and g
C = [C1; C2; C3];


% % 
% Build Vectors

% Get vectors of coefficients of f(x,y), g(x,y) and h(x,y)
fxy_vec = GetAsVector_Version1(fxy);
gxy_vec = GetAsVector_Version1(gxy);
hxy_vec = GetAsVector_Version1(hxy);

% % Build the RHS vector
rhs_vec = [fxy_vec;
           gxy_vec;
           hxy_vec];

% % Calculate the solution vector

% Calculate the x vector by pinv       
x = SolveAx_b(C,rhs_vec);
dxy_vec = x;

% Arrange dw into a matrix form based on its dimensions.
dxy = GetAsMatrix_Version1(dxy_vec,t1,t2);


end
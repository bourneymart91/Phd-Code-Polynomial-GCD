function [uxy, vxy, wxy] = GetQuotients_Relative_Bivariate_3Polys(fxy, gxy, hxy, k1, k2)
% GetQuotients(fxy,gxy,k1,k2)
%
% Given two polynomials and the knowledge of the degree of the GCD. Obtain
% the two quotient polynomials u and v.
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of the polynomial h(x,y)
%
% t1 : (Int) Degree of GCD d(x,y) with respect to x
%
% t2 : (Int) Degree of GCD d(x,y) with respect to y
%
% % Outputs
%
% uxy : (Matrix) Coefficients of polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of polynomial v(x,y)
%
% wxy : (Matrix) Coefficients of polynomial w(x,y)
%

% Get the degree of polynomial f(x,y), g(x,y) and w(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);

% Build the (k1,k2)th Sylvester matrix
Sk = BuildT_Relative_Bivariate_3Polys(fxy, gxy, hxy, k1, k2);

% % Get the index of the optimal column for removal
idx_optColumn = GetOptimalColumn_Relative(Sk);

% Having found the optimal column, obtain u and v the quotient polynomials.
Atj = Sk;
cki = Sk(:,idx_optColumn);
Atj(:,idx_optColumn) = [];

% Get the solution vector.
x_ls = SolveAx_b(Atj,cki);

% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(idx_optColumn)-1);
    -1;
    x_ls(idx_optColumn:end);
    ]  ;

% Get number of coefficients in u(x,y) and v(x,y) and w(x,y)
nCoefficients_vxy = (n1-k1+1) * (n2-k2+1);
nCoefficients_wxy = (o1-k1+1) * (o2-k2+1);

% Get the vector of coefficients of v
vec_vxy = vecx(1:nCoefficients_vxy);

vec_wxy = vecx(nCoefficients_vxy + 1 : nCoefficients_vxy + nCoefficients_wxy);
      
% Get the vector of coefficients of u
vec_uxy = (-1).*vecx(nCoefficients_vxy + nCoefficients_wxy +1 : end);
        
% % Get u(x,y) and v(x,y) and w(x,y) in matrix form
uxy = GetAsMatrix_Version1(vec_uxy, m1-k1, m2-k2);
wxy = GetAsMatrix_Version1(vec_wxy, o1-k1, o2-k2);
vxy = GetAsMatrix_Version1(vec_vxy, n1-k1, n2-k2);


end

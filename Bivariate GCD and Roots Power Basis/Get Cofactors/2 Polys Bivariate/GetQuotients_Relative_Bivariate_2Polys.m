function [uxy, vxy] = GetQuotients_Relative_Bivariate_2Polys(fxy, gxy, k1, k2)
% GetQuotients_Relative_Bivariate_2Polys(fxy, gxy, k1, k2)
%
% Given two polynomials and the knowledge of the degree of the GCD. Obtain
% the two quotient polynomials u and v.
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y) 
%
% gxy : (Matrix) Coefficients of the polynomial f(x,y) 
%
% t1 : (Int) Degree of GCD d(x,y) with respect to x
%
% t2 : (Int) Degree of GCD d(x,y) with respect to y


% Get the degree of polynomial f(x,y) and g(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

% Build the (k1,k2)th Sylvester subresultant matrix
Sk = BuildT_Relative_Bivariate_2Polys(fxy, gxy, k1, k2);

% % Get the index of the optimal column for removal from the Sylvester
% subresultant matrix
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

% Get number of coefficients in u(x,y) and v(x,y)
nCoeffs_vxy = (n1-k1+1) * (n2-k2+1);

% Get the vector of coefficients of v(x,y)
vxy_calc = vecx(1:nCoeffs_vxy);
      
% Get the vector of coefficients of u(x,y)
uxy_calc = (-1).*vecx(nCoeffs_vxy+1:end);
        

% % Get u and v in matrix form
% Arrange uw into a matrix form based on their dimensions.
uxy = GetAsMatrix(uxy_calc, m1-k1, m2-k2);
vxy = GetAsMatrix(vxy_calc, n1-k1, n2-k2);


end

function Y_kk1k2 = BuildY_BothDegree_STLN(x, m, m1, m2, n, n1, n2, k, k1, k2)
% BuildY_BothDegree_STLN(x, m, m1, m2, n, n1, n2, k, k1, k2)
%
% Build the matrix Y_{t} Where Y(x) * z = E(z) * x
% The vector x only contains the non-zero entries of v(x,y) and u(x,y)
%
% Inputs
%
% x : (Vector) obtained by inserting zero into the Least squares solution.
%
% m : (Int) Total degree of polynomial f(x,y)
%
% m1 : (Int) Degree of the polynomial f(x,y) with respect to x
%
% m2 : (Int) Degree of the polynomial f(x,y) with respect to y
%
% n : (Int) Total degree of polynomial g(x,y)
%
% n1 :(Int) Degree of the polynomial g(x,y) with respect to x
%
% n2 : (Int) Degree of the polynomial g(x,y) with respect to y
% 
% k : Total degree of polynomial d(x,y)
%
% k1 : (Int) Degree of polynomial d(x,y) with respect to x 
%
% k2 : (Int) Degree of the polynomial d(x,y) with respect to y
%
% % Outputs
%
% Y_kk1k2 : The matrix Y_{k,k_{1},k_{2}}

% Get number of coefficients in u(x,y) and v(x,y)
nCoefficients_uxy = (m1-k1+1) * (m2-k2+1);
nCoefficients_vxy = (n1-k1+1) * (n2-k2+1);

% Get number of non-zero Coefficients in u(x,y) and v(x,y)
nNonZeros_uxy = GetNumNonZeros(m1-k1,m2-k2,m-k);
nNonZeros_vxy = GetNumNonZeros(n1-k1,n2-k2,n-k);

% Get number of zero coefficients in u(x,y) and v(x,y)
nZeros_uxy = nCoefficients_uxy - nNonZeros_uxy;
nZeros_vxy = nCoefficients_vxy - nNonZeros_vxy;



% % Split the vector x
x1 = x(1:nNonZeros_vxy);
x2 = 1.* x(nNonZeros_vxy + 1 : end);

% get x_u as a matrix
mat_x2 = GetAsMatrix_Version1(...
    [
        x2;
        zeros(nZeros_uxy,1)
    ]...
    ,m1-k1,m2-k2);

% Get x_v as a matrix
mat_x1 = GetAsMatrix_Version1(...
    [
        x1;
        zeros(nZeros_vxy,1);
    ]...
    ,n1-k1,n2-k2);

% Build the matrix C1 and C2
T_x1 = BuildT1_Both_Bivariate(mat_x1, n-k, m, m1, m2);
T_x2 = BuildT1_Both_Bivariate(mat_x2, m-k, n, n1, n2);

% Build the matrix Y_{k}
Y_kk1k2 = [T_x1 T_x2];



end
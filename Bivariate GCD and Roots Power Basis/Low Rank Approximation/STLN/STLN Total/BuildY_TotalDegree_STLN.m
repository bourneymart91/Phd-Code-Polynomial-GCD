function Y = BuildY_TotalDegree_STLN(x,m,n,k)
% Build the matrix Y_{t}
% Where Y(x) * z = E(z) * x
%
% % Inputs.
%
% x : (Vector) x least squares solution with zero inserted.
%
% m : (Int) Degree of f(x,y)
% 
% n : (Int) Degree of g(x,y)
%
% k : (Int) Degree of d(x,y)
%
% % Outputs
%
% Y : (Matrix)

% Get the  number of coefficients of x_{1}
nCoefficients_x1 = nchoosek(n-k+2, 2);
nCoefficients_x2 = nchoosek(m-k+2, 2);

nZeros_x1 = nchoosek(n-k+1, 2);
nZeros_x2 = nchoosek(m-k+1, 2);

% Get vector of coefficients of x_{1}(x,y)
x1 = x(1:nCoefficients_x1);

% Get vector of coefficients of x_{2}(x,y)
x2 = x(nCoefficients_x1+1:nCoefficients_x1 + nCoefficients_x2);

% Get x_{u}(x,y) and x_{v}(x,y) as a matrix
matrix_x1 = GetAsMatrix_Version1([x1; zeros(nZeros_x1,1)], n-k, n-k);
matrix_x2 = GetAsMatrix_Version1([x2; zeros(nZeros_x2,1)], m-k, m-k);

% Build the matrices C(v) and C(u)
C1 = BuildT1_Total_Bivariate(matrix_x1, n-k, m);
C2 = BuildT1_Total_Bivariate(matrix_x2, m-k, n);

% Build the Matrix

Y = [C1 C2];



end
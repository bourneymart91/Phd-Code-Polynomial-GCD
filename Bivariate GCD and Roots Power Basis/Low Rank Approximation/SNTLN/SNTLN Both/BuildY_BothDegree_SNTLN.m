function Y = BuildY_BothDegree_SNTLN(x,m,m1,m2,n,n1,n2,k,k1,k2,alpha,th1,th2)
% This function builds the matrix Y_{t_{1},t_{2}}, 
% Where Y_{t1,t2}(x1,x2)*z = E_{t_{1},t_{2}}(z)*x
%
% % Inputs
%
% m : (Int) Total degeree of polynomial f(x,y) 
%
% m1 : (Int) Degree of f(x,y) with respect to x
%
% m2 : (Int) Degree of f(x,y) with respect to y
% 
% n : (Int) Total degree of polynomial g(x,y)
%
% n1 : (Int) Degree of g(x,y) with respect to x
%
% n2 : (Int) Degree of g(x,y) with respect to y
%
% k : (Int) Total degree of d(x,y)
%
% k1 : (Int) Degree of d(x,y) with respect to x
%
% k2 : (Int) Degree of d(x,y) with respect to y
%
% idx_col : (Int) Index of column removed from Sylvester subresultant matrix
%
% x_ls : (Vector) Least squares solution to A_{k,k1,k2}x = c_{k,k1,k2}
%
% alpha : (Float) Optimal value of \alpha
% 
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}
%
% % Outputs
%
% Y : (Matrix) Y(x1,x2)


% Separate the x into x1 and x2

nCoefficients_x1 = GetNumNonZeros(n1-k1,n2-k2,n-k);
nCoefficients_x2 = GetNumNonZeros(m1-k1,m2-k2,m-k);

% Get number of zero and non-zero coefficients of f(x,y)
nNonZeros_fxy = GetNumNonZeros(m1,m2,m);
nZeros_fxy = (m1+1)*(m2+1) - nNonZeros_fxy;

% Get number of zero and non-zero coefficients of g(x,y)
nNonZeros_gxy = GetNumNonZeros(n1,n2,n);
nZeros_gxy = (n1+1)*(n2+1) - nNonZeros_gxy;


nCoeffs_x1_mat = (n1-k1+1) * (n2-k2+1);
nCoeffs_x2_mat = (m1-k1+1) * (m2-k2+1);

nZeros_x1_mat = nCoeffs_x1_mat - nCoefficients_x1;
nZeros_x2_mat = nCoeffs_x2_mat - nCoefficients_x2;


% Get the vector of coefficients of x1
x1_vec = x(1:nCoefficients_x1);

% Get the vector of coefficients of x2
x2_vec = x(nCoefficients_x1+1:nCoefficients_x1 + nCoefficients_x2);

x1_vec = [x1_vec ; zeros(nZeros_x1_mat,1)];
x2_vec = [x2_vec ; zeros(nZeros_x2_mat,1)];

% Get x1 and x2 as matrices.
x1_mat = GetAsMatrix_Version1(x1_vec, n1-k1, n2-k2);
x2_mat = GetAsMatrix_Version1(x2_vec, m1-k1, m2-k2);

% Build the convolution matrix T_{m,m1,m2}(x1)
T1 = BuildT1_Both_Bivariate(x1_mat, n-k, m, m1, m2);

% Build the convolution matrix T_{n,n1,n2}(x2)
T2 = BuildT1_Both_Bivariate(x2_mat, m-k, n, n1, n2);

% Get matrix of ones with same structure as matrix of coefficients of
% f(x,y)
f_vec = [ones(nNonZeros_fxy,1); zeros(nZeros_fxy,1)];
f_mat = GetAsMatrix_Version1(f_vec,m1,m2);

% Get matrix of ones with same structure as matrix of coefficients of
% g(x,y)
g_vec = [ones(nNonZeros_gxy,1); zeros(nZeros_gxy,1)];
g_mat = GetAsMatrix_Version1(g_vec,n1,n2);

% Get a vector of theta_{1}theta_{2} corresponding to entries of the
% polynomial f(x,y)
th_f = GetWithThetas(f_mat,th1,th2);
th_f = GetAsVector_Version1(th_f);
th_f = th_f(1:nNonZeros_fxy);
th_f = diag(th_f);

% Get a vector of \theta_{1}\theta_{2} corresponding to entires of the
% polynomial g(x,y)
th_g = GetWithThetas(g_mat,th1,th2);
th_g = GetAsVector_Version1(th_g);
th_g = th_g(1:nNonZeros_gxy);
th_g = diag(th_g);


% multiply by the alpha of g
Y = [T1*th_f alpha.*T2*th_g] ;
end
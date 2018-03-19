function Y = BuildY_RelativeDegree_SNTLN(x,m1,m2,n1,n2,k1,k2,alpha,th1,th2)
% This function builds the matrix Y_{t_{1},t_{2}}, 
% Where Y(x)*z = E_{t_{1},t_{2}}(z)*x
%
% % Inputs
%
% x : 
%
% m1 : Degree of polynomial f(x,y) with respect to x
%
% m2 : Degree of polynomial f(x,y) with respect to y
%
% n1 : Degree of polynomial g(x,y) with respect to x
%
% n2 : Degree of polynomial g(x,y) with respect to y
%
% k1 : Degree of polynomial d(x,y) with respect to x
%
% k2 : Degree of polynomial d(x,y) with respect to y
%
% alpha : Optimal value of \alpha
%
% th1 : Optimal value of \theta_{1}
% 
% th2 : Optimal value of \theta_{2}
%
% % Outputs
%
% Y : Matrix Y


% Separate the x into x1 and x2.

% Get number of coefficients in x_{1}
nCoeffs_x1 = (n1-k1+1) * (n2-k2+1);

% Get vector of coefficients of x_{1}
x1_vec = x(1:nCoeffs_x1);

% Get the vector of coefficients of x_{2}
x2_vec = x(nCoeffs_x1+1:end);

% Get the matrix of coefficients of x_{1}
x1_xy = GetAsMatrix(x1_vec,n1-k1,n2-k2);

% Get the matrix of coefficients of x_{2}
x2_xy = GetAsMatrix(x2_vec,m1-k1,m2-k2);

% Build convolution matrix of x_{1}
T1 = BuildT1_Relative_Bivariate(x1_xy, m1, m2);

% Build convolution matrix of x_{2}
T2 = BuildT1_Relative_Bivariate(x2_xy, n1, n2);

% Get a vector of thetas corresponding to coefficients of f(\omega,\omega)
th1_mat = diag(th1.^(0:1:m1));
th2_mat = diag(th2.^(0:1:m2));
fww_thetas_mat = th1_mat * ones(m1+1,m2+1) * th2_mat;
th_fww = GetAsVector(fww_thetas_mat);
th_fww = diag(th_fww);

% Get a vector of thetas corresponding to coefficients of g(\omega,\omega)
th1_mat = diag(th1.^(0:1:n1));
th2_mat = diag(th2.^(0:1:n2));
gww_thetas_mat = th1_mat * ones(n1+1,n2+1) * th2_mat;
th_gww = GetAsVector(gww_thetas_mat);
th_gww = diag(th_gww);

% Build the matrix Y
Y = [T1*th_fww alpha*T2*th_gww];
end
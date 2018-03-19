function DYQ = BuildDYQ_SNTLN(xk, m, n, k, alpha, theta)
% USED IN SNTLN function
%
% Construct Matrix DYQ, such that E_{k}(z)x = D^{-1}Y_{k}(x)Qz, where E_{k}(z) is a
% matrix of structured perturbations applied to S_{k}, where S_{k} = DTQ.
%
% Inputs
%
% xk : (Vector)
%
% m : (Int) Degree of polynomial f(x)
%
% n : (Int) Degree of polynomial g(x)
%
% k : (Int) Degree of GCD
%
% alpha : (Float) The optimal value of \alpha
%
% theta : (Float) The optimal value of \theta
%
% % Outputs
%
% DYQ : (Matrix) 

% Get the number of coefficients in x1
nCoefficients_x1 = n - k + 1;

% Split the vector x_{k} = [x_{1} x_{2}]
x1 = xk(1 : nCoefficients_x1);
x2 = xk(nCoefficients_x1 + 1 : end);

% Build the matrix D^{-1}_{m+n-k}
D = BuildD_2Polys(m, n - k);

% Build the matrices Y_{1} and Y_{2}
Y1 = BuildT1(x1, m);
Y2 = BuildT1(x2, n);

% Build the matrices Q_{m} and Q_{n}
Qm = BuildQ1(m);
Qn = BuildQ1(n);

% Get a diagonal matrices of thetas corresponding to polynomials f(x) and
% g(x).
th_f = diag(theta.^(0 : 1 : m));
th_g = diag(theta.^(0 : 1 : n));


DYQ = D*[Y1*Qm*th_f  alpha.*Y2*Qn*th_g] ;
   




end








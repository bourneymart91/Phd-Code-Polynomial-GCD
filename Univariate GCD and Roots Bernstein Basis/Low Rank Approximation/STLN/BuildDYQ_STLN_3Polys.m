function DYQ = BuildDYQ_STLN_3Polys(xk, m, n, o, k)
% USED IN SNTLN function
%
% Construct Matrix DYQ, such that E_{k}(z)x = D^{-1}Y_{k}(x)Qz, where E_{k}(z) is a
% matrix of structured perturbations applied to S_{k}, where S_{k} = DTQ.
%
% Inputs
%
% xk : (Vector) Vector
%
% m : (Int) Degree of polynomial f(x,y)
%
% n : (Int) Degree of polynomial g(x,y)
%
% k : (Int) Degree of GCD
%
% alpha : (Float) The optimal value of \alpha
%
% theta : (Float) The optimal value of \theta

% Get the number of coefficients in x1
nCoefficients_x1 = n - k + 1;
nCoefficients_x2 = o - k + 1;

% Split the vector x_{k} = [x_{1} x_{2}]
x1 = xk(1 : nCoefficients_x1);
x2 = xk(nCoefficients_x1 + 1 : nCoefficients_x1 + nCoefficients_x2);
x3 = xk(nCoefficients_x1 + nCoefficients_x2 + 1 : end);


% Build the matrix D^{-1}_{m+n-k}
D = BuildD_3Polys(m, n - k, o - k);

% Build the matrices Y_{1} and Y_{2}
Y1 = BuildT1(x1, m);
Y2 = BuildT1(x3, n);
Y3 = zeros(m + n - k + 1, o + 1);
Y4 = BuildT1(x2, m);
Y5 = zeros(m + o - k + 1, n + 1);
Y6 = BuildT1(x3, o);


% Build the matrices Q_{m} and Q_{n}
Qm = BuildQ1(m);
Qn = BuildQ1(n);
Qo = BuildQ1(o);
Q = blkdiag(Qm, Qn, Qo);

% Build DYQ
DYQ = D*[Y1  Y2 Y3 ; Y4 Y5 Y6] * Q ;


end

function DYQ = BuildDYQ_STLN_RootSpecific(xk, m, n, k, alpha, theta)
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
nCoefficients_x1 = m - k;

% Split the vector x_{k} = [x_{1} x_{2}]
x1 = xk(1 : nCoefficients_x1);
x2 = xk(nCoefficients_x1 + 1 : end);

% Build First Partition

D = BuildD_2Polys(m, m - k - 1);

th_m = diag(GetWithThetas(ones(m+1,1),theta));
Tx1 = BuildT1(x1, m);
Qm = BuildQ1(m);

FirstPart = D* Tx1 * th_m * Qm;

% Build Second Partition


Tx2 = BuildT1(x2, m - 1);

D = BuildD_2Polys(m, m - k - 1);

Qn = BuildQ1(m-1);

th_n = diag(GetWithThetas(ones(m,1),theta));

mat1 = ...
    [
        zeros(2*m - k ,1) alpha.* Tx2 * Qn * th_n
    ];

mat2 = ...
    [
        alpha .* Tx2 * Qn * th_n zeros(2*m - k, 1)
    ];

SecondPart = D * m*(mat1 - mat2);


DYQ = [FirstPart + SecondPart];


end




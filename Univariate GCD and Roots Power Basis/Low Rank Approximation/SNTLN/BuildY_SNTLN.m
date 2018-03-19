function Y = BuildY_SNTLN(m, n, k, x, alpha, theta)
% Build the Matrix Y such that Y(x)*z = E(z)*x
%
% % Inputs
%
% m : (Int) Degree of polynomial f(x)
%
% n : (Int) Degree of polynomial g(x)
%
% k : (Int) Degree of polynomial d(x)
%
% x : (Vector) x = [x1 x2]
%
% alpha : (Float) Optimal Value of \alpha
%
% theta : (Float) Optimal value of \theta

% Get the vectors u(\omega,\theta) and v(\omega,\theta)

% The matrix Y is such that S_{k}(f,g)*x = Y(x)*[f;g]
nCoefficients_x1 = n-k+1;

% Get the x values corresponding to v_{k}
x1 = x(1:nCoefficients_x1);

% Get the x values corresponding to u_{k}
x2 = x(nCoefficients_x1 + 1:end);


T1_x1 = BuildT1(x1, m);
T2_x2 = BuildT1(x2, n);

% Get thetas corresponding to polynomial f(x)
th_f = diag(GetWithThetas(ones(m+1, 1), theta));

% Get thetas corresponding to polynomial g(x)
th_g = diag(GetWithThetas(ones(n+1, 1), theta));

% Build the matrix Y
Y = [T1_x1*th_f alpha.*T2_x2*th_g];


end

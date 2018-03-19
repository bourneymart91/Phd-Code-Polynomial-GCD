function Q = BuildQ_3Polys(m, n, o, k)
% Build the diagonal matrix Q corresponding to the binomial coefficients
% of coprime polynomials u, v and w.
%
%
% Inputs
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% o : Degree of polynomial h(x)
%
% k : Index of subresultant S_{k} to be formed.
%
%
% Outputs.
%
% Q : The diagonal matrix of binomial coefficients corresponding to coprime
%       polynomials u and v.
%
%

% Build first partition of Q corresponding to the binomial coefficients of
% v(y). \binom{n-k}{i}
Q1 = BuildQ1(n - k);

% Build second partition of Q
Q2 = BuildQ1(o - k);

% Build third partition of Q corresponding to the binomial coefficients of
% u(y). \binom{m-k}{i}
Q3 = BuildQ1(m - k);


% Join the two partitions as a diagonal matrix.
Q = blkdiag(Q1, Q2, Q3);



end

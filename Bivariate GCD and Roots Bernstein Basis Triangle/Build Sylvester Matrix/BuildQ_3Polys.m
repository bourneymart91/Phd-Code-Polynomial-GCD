function Q = BuildQ_3Polys(m, n, o, k)
% BuildQ_3Polys(m, n, o, k)
%
% Build the block diagonal matrix Q_{k}, part of the modified Sylvester
% subresultant matrix S_{k}(f,g). Q_{k} is the block diagonal matrix
% consisting of three submatrices Q_{n-k} Q_{o-k} and Q_{m-k}
%
% Inputs
%
% [m, n, o] : Degree of polynomial f(x,y), g(x,y) and h(x,y)
%
% k : Index of subresultant matrix being constructed.

% Build the matrix Q_{n-k} 
Q1 = BuildQ1(n-k);

% Build the matrix Q_{o-k}
Q2 = BuildQ1(o-k);

% Build the matrix Q_{m-k}
Q3 = BuildQ1(m-k);

% Build the matrix Q_{k}
Q = blkdiag(Q1,Q2,Q3);

end
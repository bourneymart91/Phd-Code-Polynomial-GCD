function Q = BuildQ_2Polys(m, n, k)
% BuildQ_2Polys(m, n, k)
%
% Build the block diagonal matrix Q_{k}, part of the modified Sylvester
% subresultant matrix S_{k}(f,g). Q_{k} is the block diagonal matrix
% consisting of two submatrices Q_{n-k} and Q_{m-k}
%
% Inputs
%
% m : Degree of polynomial f(x,y)
%
% n : Degree of polynomial g(x,y)
%
% k : Index of subresultant matrix being constructed.

% Build the matrix Q_{n-k} 
Q1 = BuildQ1(n-k);

% Build the matrix Q_{m-k}
Q2 = BuildQ1(m-k);

% Build the matrix Q_{k}
Q = blkdiag(Q1,Q2);

end
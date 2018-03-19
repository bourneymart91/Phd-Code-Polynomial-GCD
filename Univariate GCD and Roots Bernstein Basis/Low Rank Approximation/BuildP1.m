function P1 = BuildP1(m, n_k, idx_col)
% Build the matrix P1, a partition of P, where P*[f;g] gives a column of
% the Sylvester subresultant matrix S_{k}(f,g)
%
% % Inputs
%
% m : (Int) Degree of polynomial f(x)
%
% n_k : (Int) Degree of polynomial v(x)
%
% idx_col : (Int) Index of column of T_{n-k}(f) which is being constructed.
%
% % Outputs
%
% P1 : (Matrix) Partition of the matrix P

% Initialise the partition P_{1}
P1 = zeros(m+n_k+1,m+1);

% Get with binomials corresponding to v(x)
mat = eye(m+1);

% 
P1(idx_col : idx_col + m, :) = mat;

end

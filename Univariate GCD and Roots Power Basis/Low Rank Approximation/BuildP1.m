function P1 = BuildP1(m,n,k,idx_col)
% BuildP1_SNTLN(m,n,k,idx_col)
% 
% Build the matrix P1
% 
% % Inputs
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% k : Degree of divisor d(x) and index of kth subresultant S_{k}(f,g)
%
% % Outputs
% 
% P1 : Partition of the matrix P.

% P1 is a diagonalisation of a vector given by [zeros; ones; zeros]
num_zero_rows_top = idx_col-1;
num_zero_rows_bottom = (m+n-k+1) - (m+1) - num_zero_rows_top;
P1 = ...
    [
    zeros(num_zero_rows_top,m+1);
    eye(m+1);
    zeros(num_zero_rows_bottom,m+1);
    ];
end
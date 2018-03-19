function [W] = BuildW(m,n,k)
% Build the matrix W.
%
% % Inputs
% 
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% k : Degree of polynomial d(x)

W1 = BuildWPartition(n-k);
W2 = BuildWPartition(m-k);

W = blkdiag(W1,W2);

end

function [W1] = BuildWPartition(n_k)

vec = (n_k+1 : -1 : 1);

vec = vec ./ (n_k+1);

W1 = [...
        diag(vec) ; 
        zeros(1,n_k+1)
     ];

end
function D = BuildD_3Polys_2Eqns(m, n_k, o_k)
% Build the matrix D^{-1} which is the diagonal matrix of binomial 
% coefficients. Used to construct the Sylvester Subresultant matrices in 
% format D_{k}^{-1}*T_{k}(f,g)*Q.
% 
%
%
% % Input
%
% m : Degree of polynomial f(x)
%
% n_k : Degree of polynomial v(x)
%
% o_k : Degree of polynomial w(x)
%
% % Outputs
%
% D : The matrix D^{-1} = [D^{-1}_{m+n-k}       0,        ]
%                         [     0           D^{-1}_{m+o-k}]

D1 = BuildD_2Polys(m, n_k);
D2 = BuildD_2Polys(m, o_k);

D = blkdiag(D1, D2);

end
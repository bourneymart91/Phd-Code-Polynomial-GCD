function D = BuildD_3Polys_3Eqns(m, n, o, k)
% Build the matrix D^{-1} which is the diagonal matrix of binomial 
% coefficients. Used to construct the Sylvester Subresultant matrices in 
% format D_{k}^{-1}*T_{k}(f,g)*Q.
% 
%
%
% % Input
%
% m : (Int) Degree of polynomial f(x)
%
% n : (Int) Degree of polynomial g(x)
%
% o : (Int) Degree of polynomial h(x)
%
% k : (Int)
%
% % Outputs
%
% D :

D1 = BuildD_2Polys(m, n - k);
D2 = BuildD_2Polys(m, o - k);
D3 = BuildD_2Polys(o, n - k);

D = blkdiag(D1, D2, D3);

end
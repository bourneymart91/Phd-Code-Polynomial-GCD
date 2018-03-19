function D = BuildD_3Polys_3Eqns(m, n, o, k)
% BuildD(m,n_t)
%
% Build the diagonal matrix D^{-1} for the convolution of two polynomials
% f(x,y) and v(x,y) of degrees m and n-t.
%
% Inputs
%
% m : (Int) Degree of polynomial f(x,y)
%
% n : (Int) Degree of polynomial g(x,y)
%
% o : (Int)
%
% k : (Int)

D1 = BuildD_2Polys(m, n - k);

D2 = BuildD_2Polys(m, o - k);

D3 = BuildD_2Polys(n, o - k);

D = blkdiag(D1, D2, D3);

end

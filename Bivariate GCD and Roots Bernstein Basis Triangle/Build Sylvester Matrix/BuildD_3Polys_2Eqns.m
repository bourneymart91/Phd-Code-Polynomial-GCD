function D = BuildD_3Polys_2Eqns(m,n,o,k)
% BuildD(m,n_t)
%
% Build the diagonal matrix D^{-1} for the convolution of two polynomials
% f(x,y) and v(x,y) of degrees m and n-t.
%
% Inputs
%
% m : (Int) Degree of polynomial f(x,y)
%
% n_t : (Int) Degree of polynomial v(x,y)

D1 = BuildD_2Polys(m, n-k);

D2 = BuildD_2Polys(m, o-k);

D = blkdiag(D1, D2);

end

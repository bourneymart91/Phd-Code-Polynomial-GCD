function [] = o_GCD_matlab(fx_brn,gx_brn)

% Given 
% Get the Binomial coefficients corresponding to the coefficients of
% polynomial f.
fx_bi = GetWithBinomials(fx_brn);
gx_bi = GetWithBinomials(gx_brn);

% % Matlab roots function expects descending powers of x^n ... x 1
fx_bi = flipud(fx_bi);
gx_bi = flipud(gx_bi);


[d,u,v] = gcd(fx_bi, gx_bi)



end
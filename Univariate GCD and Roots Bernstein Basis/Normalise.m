function fx_n = Normalise(fx)
% Given the column vector of coefficients of polynomial f(x). Normalise by
% dividing by the first coefficient.

fx_n = fx./ fx(1);

end
function fx = GetWithoutThetas(fw,theta)
% Given the coefficients of the polynomial f(w), get the coefficients of
% f(x), where x = \theta\omega

% Get the degree of polynomial f(x)
m = GetDegree(fw);

% Get coefficients of f(x) from f(w)
fx = fw./(theta.^(0:1:m))';

end
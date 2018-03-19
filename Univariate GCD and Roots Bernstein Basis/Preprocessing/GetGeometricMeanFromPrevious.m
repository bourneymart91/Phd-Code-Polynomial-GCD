function [GM] = GetGeometricMeanFromPrevious(fx, GM_prev, m, n_k)
% This function allows for the computation of the geometric mean of the
% entries in the $k$th subresultant matrix when given the geometric mean of
% the entries in the (k-1)th subresultant matrix.
%
% % Inputs
%
% fx : (Vector) The coefficients of the polynomial f(x)
%
% GM_Prev : (float) The Geometric mean of coefficients of f(x) in the 
% previous Sylvester subresultant matrix S_{k-1}
%
% m : (Int) The degree of polynomial f(x)
%
% n_k : (Int) is the degree of polynomial v(x) and determines the number of columns in
% T_{n-k}(f(x)), where nCols = n-k+1.


p1 = (m+n_k+1) ./ (n_k+1);

p1 = p1 .^((n_k+2)./(n_k+1));

p2 = GM_prev .^((n_k+2)./(n_k+1));

p3 = nchoosek(m+n_k,m) .^ (1./(n_k+1));

prod_a = 1;
prod_b = 1;

for i = 0:1:m

    prod_a = prod_a * abs(fx(i+1));
    prod_b = prod_b * nchoosek(n_k+1+i,i);

end

p4 = (1./(prod_a * prod_b * prod_b)) .^(1./((n_k+1)*(m+1)));

GM = p1 * p2 * p3 * p4;

end

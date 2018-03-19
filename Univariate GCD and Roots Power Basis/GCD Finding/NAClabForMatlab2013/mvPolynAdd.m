% Compute the addition f+g of two multivariate polynomials f and g, 
%    or the linear combination  s*f + t*g  if scalars s and t are also given
%
%    Syntax:  >> d = mvPolynAdd(f,g)
%             >> d = mvPolynAdd(f,g,s,t)
%
%    Input:   f, g -- multivariate polynomials in coeff. matrices
%             s, t -- (optional) complex numbers
%
%   Output:   multivariate polynomial  f+g,  or s*f + t*g

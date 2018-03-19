function fw = GetWithThetas(fx, th)
% Take the coefficients of polynomial f(x,y) and include thetas so that
% the ouput is the coefficients of f(w,w), in the modified Bernstein
% basis.

m = GetDegree(fx);

% Multiply coefficients of f(x) a_{i} by \theta^{i}
fw = fx .* (th.^(0:1:m)');


end
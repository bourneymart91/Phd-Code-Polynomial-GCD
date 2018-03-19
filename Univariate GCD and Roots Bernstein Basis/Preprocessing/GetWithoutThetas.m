function fw = GetWithoutThetas(fx,th)
% Take the coefficients of polynomial f(\omega,\theta), divide each 
% coefficient a_{i} by \theta^{i} to obtain f(x)


% Get the degree of the polynomial f(x)
m = GetDegree(fx);

% Divide coefficients of f(x),  a_{i} by \theta^{i}
fw = fx ./ (th.^(0:1:m)');


end
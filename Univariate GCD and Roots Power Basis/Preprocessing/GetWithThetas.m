function fw = GetWithThetas(fx,th)
% Given coefficients of the polynomial f(x) in the power basis, multiply
% the coefficient a_{i} by \theta^{i} to obtain coefficients in the scaled
% basis. where x = \theta \omega.
%
% % Inputs
%
% fx : Coefficients of the polynomial f(x)
%
% th : Optimal value of \theta


% Get the degree of polynomial f(x)
m = GetDegree(fx);

% Calculate coefficients of f(w).
fw = diag(th.^(0:1:m)) * fx;

end
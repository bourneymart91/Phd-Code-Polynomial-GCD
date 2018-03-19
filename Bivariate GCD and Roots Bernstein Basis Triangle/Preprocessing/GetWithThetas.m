function [fww] = GetWithThetas(fxy,m,th1,th2)
% Get f(\omega_{1},\omega_{2})
%
% % Inputs
%
% fxy : Coefficients of polynomial f(x,y)
%
% m : Total degree of polynomial f(x,y)
%
% th1 : \theta_{1}
%
% th2 : \theta_{2}


% Multiply the rows of fxy by th_{1}^{i}
pre_thetas = diag(th1.^(0:1:m));

% Multiply the columns of the matrix by th_{2}^{j}
post_thetas = diag(th2.^(0:1:m));

% Get f(\omega_{1},\omega_{2})
fww = pre_thetas * fxy * post_thetas;



end
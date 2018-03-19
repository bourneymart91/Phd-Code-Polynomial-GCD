function fxy = GetWithoutThetas(fww, th1, th2)
% Given the coefficients of the polynomial f(w,w), divide the (i,j)-th 
% entry by theta_{1}^{i}\theta_{2}^{j}.
%
% % Inputs
%
% fww : Coefficients of polynomial f(\omega_{1},\omega_{2})
%
% th1 : \theta_{1}
%
% th2 : \theta_{2}
%
% % Outputs
% 
% fxy : Coefficients of polynomial f(x,y)

% Get the degree of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fww);

% Calculate f(x,y) from f(w,w)

% Build matrix th1_mat
th1_mat = diag(1./(th1.^(0:1:m1)));

% Build matrix th2_mat
th2_mat = diag(1./(th2.^(0:1:m2)));

% Multiply coefficient matrix of f(\omega_{1},\omega_{2}) by the theta
% matrices to get coefficient matrix of f(x,y)
fxy = th1_mat * fww * th2_mat;

end
function fww = GetWithThetas(fxy,th1,th2)
% GetWithThetas(fxy,th1,th2)
%
% Inputs.
%
% fxy : Matrix of coefficients of polynomial f(x,y)
%
% th1 : \theta_{1}
%
% th2 : \theta_{2}
%
% % Outputs
%
% fww : Matrix of coefficients of polynomial f(\omega_{1},\omega_{2})

% Get degree of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Multiply row i by theta_{1}.^{i} and columns j by \theta_{2}.^{j}
fww = diag(th1.^(0:1:m1)) * fxy * diag(th2.^(0:1:m2));


end
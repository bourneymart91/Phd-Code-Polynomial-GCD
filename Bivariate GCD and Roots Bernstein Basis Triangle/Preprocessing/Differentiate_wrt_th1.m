function fww_wrt_th1 = Differentiate_wrt_th1(fww,th1)
% Get the partial derivative of the polynomial f(\omega_{1},\omega_{2}) 
% with respect to \theta_{1}
%
% % Inputs.
%
% fww : (Matrix) Coefficients of f(\omega_{1},\omega_{2})
%
% th1 : (Float) \theta_{1}
%
% % Outputs
%
% fww_wrt_th2 : (Matrix) Coefficients of the partial derivative of
%               f(\omega_{1},\omega_{2}) with respect to \theta_{2}

% Get the degree of f(w,w) with respect to x
[m1,~] = GetDegree_Bivariate(fww);

temp_mat = diag((0:1:m1)./th1);

fww_wrt_th1 = temp_mat * fww;

end
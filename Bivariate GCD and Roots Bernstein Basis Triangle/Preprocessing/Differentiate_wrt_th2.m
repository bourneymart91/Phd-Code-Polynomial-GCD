function fww_wrt_th2 = Differentiate_wrt_th2(fww,th2)
% Get the partial derivative of the polynomial f(\omega_{1},\omega_{2}) 
% with respect to \theta_{2}
%
% % Inputs.
%
% fww : (Matrix) Coefficients of f(\omega_{1},\omega_{2})
%
% th2 : (Float) \theta_{2}
%
% % Outputs
%
% fww_wrt_th2 : (Matrix) Coefficients of the partial derivative of
%               f(\omega_{1},\omega_{2}) with respect to \theta_{2}


% Get the degree of f(w,w) with respect to y.
[~,m2] = GetDegree_Bivariate(fww);

temp_mat = diag((0:1:m2)./th2);

fww_wrt_th2 =  fww * temp_mat;
    
end
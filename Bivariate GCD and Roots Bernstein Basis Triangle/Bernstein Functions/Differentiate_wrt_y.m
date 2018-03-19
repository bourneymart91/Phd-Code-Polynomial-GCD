function [d_fxy] = Differentiate_wrt_y(fxy)
% Differentiate_wrt_y()
%
% Differentiate the polynomial f(x,y) with respect to x.
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% % Outputs.
%
% d_fxy : (Matrix) The partial derivative of f(x,y) with respect to x


fxy = transpose(fxy);

d_fxy = Differentiate_wrt_x(fxy);

d_fxy = transpose(d_fxy);


end
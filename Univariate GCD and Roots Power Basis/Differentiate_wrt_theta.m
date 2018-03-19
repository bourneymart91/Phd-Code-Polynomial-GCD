function f_wrt_theta = Differentiate_wrt_theta(fw,theta)
% Given the polynomial f(\omega,\theta) differentiate with respect to theta
%
% % Inputs
%
% fw : (Vector) Coefficients of polynomial f(\omega)
%
% % Outputs
%
% f_wrt_theta : (Vector) Coefficients of derivative of f(\omega) with
% respect to theta.

% Get the degree of polynomial f(\theta)
m = GetDegree(fw);

% Get transformation matrix
mat = diag(0:1:m) ./ theta;

% Get derivative of f(\theta) with respect to theta 
f_wrt_theta = mat * fw;


end
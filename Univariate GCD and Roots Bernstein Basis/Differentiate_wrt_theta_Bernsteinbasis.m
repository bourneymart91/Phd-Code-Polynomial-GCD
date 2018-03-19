function Partial_fw_wrt_theta = Differentiate_wrt_theta_Bernsteinbasis(fw, theta)
% Differentiate f(\omega,\theta) with respect to theta
%
% % Inputs
%
% fw : (Vector) Coefficients of the polynomial f(\omega)
%
% theta : (Float) Optimal value of \theta
%
% % Outputs
%
% Partial_fw_wrt_theta : (Vector) Partial derivative of f(\omega) with
% respect to \theta.



% Get the degree of polynomial f(\omega)
m = GetDegree(fw);

% Get vector [0,a0,a1,...a_{m-1}]
vec1 = fw(1:m+1);

% Get vector [0,a2,a3,...,a_{m}]
vec2 = [0; fw(1:m)];

vec3 = (vec1 - vec2).*(0:1:m)';

Partial_fw_wrt_theta = vec3./theta;
end
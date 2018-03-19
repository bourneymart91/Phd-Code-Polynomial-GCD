function [GM_fx, GM_gx, GM_hx, lambda, mu, rho, theta] = ...
    Preprocess_3Polys_2Eqns(fx, gx, hx, k)
% Preprocess_3Polys(fx, gx, hx, k)
%
% Get the optimal values of lamdba, mu, alpha and theta.
%
% % Inputs.
%
%
% fx : (Vector) Vector of coefficients of f(x)
%
% gx : (Vector) Vector of coefficients of g(x)
%
% hx : (Vector) Vector of coefficients of h(x)
%
% k : (Int) Degree of common divisor d_{k}(x)
%
% % Outputs.
%
% GM_fx : (Float) Geometric mean of entries of f(x) in k-th subresultant matrix
%
% GM_gx : (Float) Geometric mean of entries of g(x) in k-th subresultant matrix
%
% GM_hx : (Float) Geometric mean of entries of h(x) in k-th subresultant matrix
%
% alpha : (Float) Optimal value of \alpha
%
% beta : (Float) Optimal value of \beta
%
% theta : (Float) Optimal value of \theta

% Global variables
global SETTINGS

% Check number of input arguments
if (nargin ~= 4)
    error('Not enough input arguments');
end

% Get degree of polynomials f(x), g(x) and h(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

% Get the mean of the entries of f(x) in T_{n-k}(f) and T_{o-k}(f)
GM_fx = GetMean_3Polys_3Eqns(fx, n - k, o - k);
GM_gx = GetMean(gx, m - k);
GM_hx = GetMean(hx, m - k);

% Normalize f(x) and g(x) by geometric means
fx_n = fx ./ GM_fx;
gx_n = gx ./ GM_gx;
hx_n = hx ./ GM_hx;


if (SETTINGS.BOOL_ALPHA_THETA)
    
    % For each coefficient ai of F, obtain the max and min such that F_max =
    % [max a0, max a1,...] and similarly for F_min, G_max, G_min
    
    % Get maximum and minimum entries of each a_{i} in the first
    % partition T_{n-k}(f)
    [v_F_max1, v_F_min1] = GetMaxMin(fx_n, n - k);
    [v_F_max2, v_F_min2] = GetMaxMin(fx_n, o - k);
    [v_G_max, v_G_min] = GetMaxMin(gx_n, m - k);
    [v_H_max, v_H_min] = GetMaxMin(hx_n, m - k);
    
    
    % Calculate the optimal value of alpha and theta for the kth
    % subresultant matrix.
    [lambda, mu, rho, theta] = ...
        OptimalAlphaBetaGammaTheta_3Polys_2Eqns(v_F_max1, v_F_min1, ...
        v_F_max2, v_F_min2, v_G_max, v_G_min, v_H_max, v_H_min);
    
    
    
else
    lambda = 1;
    mu = 1;
    rho = 1;
    theta = 1;
    
end

end

function [GM_fx, GM_gx, GM_hx, lambda, mu, rho, theta] = Preprocess_3Polys_3Eqns(fx, gx, hx)
% Preprocess(fx, gx, hx)
%
% Preprocess the input polynomials to obtain geometric means of the
% coefficients of f(x) and g(x), and optimal values of alpha and theta such
% that the ratio of entries in sylvester matrix S(f,g) is minimized.
%
% Inputs.
%
% fx : (Vector) Vector of coefficients of polynomial f(x)
%
% gx : (Vector) Vector of coefficients of polynomial g(x)
%
% hx : (Vector) Vector of coefficients of polynomial h(x)
%
% % Outputs
%
% [GM_fx, GM_gx, GM_hx] : [Float Float Float] Geometric mean of entries 
% of f(x), g(x) and h(x)
%
% alpha : (Float)
% 
% beta : (Float)
%
% theta : (Float)

global SETTINGS

% Get Mean of entries of f(x), g(x) and h(x)
GM_fx = GetMean_2Partitions(fx);
GM_gx = GetMean_2Partitions(gx);
GM_hx = GetMean_2Partitions(hx);


% Divide f(x) and g(x) by geometric mean
fx_n = fx./ GM_fx;
gx_n = gx./ GM_gx;
hx_n = hx./ GM_hx;

switch SETTINGS.BOOL_ALPHA_THETA
    case true
        
        % Get opitmal values of alpha, beta and theta        
        [lambda, mu, rho, theta] = ...
            GetOptimalAlphaBetaAndTheta_3Polys(fx_n, gx_n, hx_n);
        

    case false
        
        % Dont apply preprocessing
        theta = 1;
        lambda = 1;
        mu = 1;
        rho = 1;
        
    otherwise
        
        error('bool_preproc either y or n');
end

end



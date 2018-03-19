function [fx_o, gx_o, dx_o, ux_o, vx_o, alpha_o, theta_o, t , GM_fx, GM_gx, rank_range] =...
    o_gcd_mymethod_Univariate_2Polys(fx, gx, limits_t, rank_range)
% o_gcd_mymethod(fx, gx, limits_t, rank_range)
%
% Given two polynomials f(x) and g(x) return the GCD d(x)
%
% % Inputs.
%
% fx : (Vector) Coefficients of polynomial f(x).
%
% gx : (Vector) Coefficients of polynomial g(x).
%
% limits_t : [(Int) (Int)] The interval which contains t
%
% % Outputs
% 
% fx_o : (Vector) Coefficients of polynomial f(x)
%
% gx_o : (Vector) Coefficients of polynomial g(x)
%
% dx_o : (Vector) Coefficients of polynomial d(x)
%
% ux_o : (Vector) Coefficients of polynomial u(x)
% 
% vx_o : (Vector) Coefficients of polynomial v(x)
%
% alpha_o : (Float) Optimal value of \alpha
%
% theta_o : (Float) Optimal value of \theta
%
% rank_range : [Float Float]



% Preprocess the polynomials f(x) and g(x)
[GM_fx, GM_gx, alpha, theta] = Preprocess_2Polys(fx, gx);

% Get f(x) and g(x) normalised by corresponding means
fx_n = fx./ GM_fx;
gx_n = gx./ GM_gx;

% Get f(\omega) from f(x) and g(\omega) from g(x)
fw = GetWithThetas(fx_n, theta);
gw = GetWithThetas(gx_n, theta);

% Get the degree of the GCD with limits defined
[t, rank_range] = GetGCDDegree_2Polys(fw, alpha.*gw, limits_t, rank_range);

% Print the degree of the GCD
LineBreakLarge();
fprintf([mfilename ' : ' sprintf('Degree of GCD : %s \n',int2str(t))]);

if (t == 0)
    dx_o = 1;
    ux_o = fx;
    vx_o = gx;
    fx_o = fx;
    gx_o = gx;
    
    alpha_o = alpha;
    theta_o = theta;
    
    return
end


% %
% %
% Get the low rank approximation and refined values for fx, gx, alpha, and
% theta. Note that vectors of polynomial coefficients are suffixed with 
% 'lr' (lr = low rank) and these are the coefficients after low rank
% approximation is computed.

[fx_lr, gx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = LowRankApprox_2Polys(fx_n, gx_n, alpha, theta, t);


[ux_lra, vx_lra, fx_lra, gx_lra, dx_lra, alpha_lra, theta_lra] = ...
    APF_2Polys(ux_lr, vx_lr, fx_lr, gx_lr, alpha_lr, theta_lr, t);

% Get f(x) and g(x)
fx_o = fx_lra;
gx_o = gx_lra;

% Get u(x) and v(x)
ux_o = ux_lra;
vx_o = vx_lra;

% Get d(x)
dx_o = dx_lra;

% Get \alpha and \theta
alpha_o = alpha_lra;
theta_o = theta_lra;
end





function [fx_o, gx_o, hx_o, dx_o, ux_o, vx_o, wx_o, ...
    lambda_o, mu_o, rho_o, theta_o, t ] = ...
    o_gcd_3Polys_mymethod(fx, gx, hx, limits_t, rank_range)
% This function computes the GCD d(x) of two noisy polynomials f(x) and g(x).
%
% Inputs:
%
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x)
%
% hx : (Vector) Coefficients of the polynomial h(x)
%
% limits_t : [Int Int] Upper and lower limits for GCD degree may be defined here
% otherwise set to [0, min(m,n)] if this is a stand alone GCD problem.
%
% rank_range : [Float Float] :  
%
%
% % Outputs:
%
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x)
%
% hx : (Vector) Coefficients of the polynomial h(x)
%
% dx : (Vector) The GCD of f(x) + \delta f(x) and g(x) + \delta g(x)
%
% ux : (Vector) Coefficients of the polynomial u(x) = f(x)/d(x)
%
% vx : (Vector) Coefficients of the polynomial v(x) = g(x)/d(x)
%
% wx : (Vector) Coefficients of the polynomial w(x) = h(x)/d(x)
%
% lambda : (Float) The optimal value \lambda_{t} obtained by preprocessing 
% the t-th subresultant matrix
%
% mu : (Float) The optimal value \mu_{t} obtained by preprocessing the t-th
% subresutlant matrix.
%
% rho : (Float) The optimal value \rho_{t} obtained by preprocessing the 
% t-th subresultant matrix.
%
% theta : (Float) The optimal value \theta_{t} obtained by preprocessing 
% the t-th subresultant matrix.






% % Get the degree of the GCD of f(x), g(x) and h(x) and values
% \lambda_{t}, \mu_{t}, \rho_{t}, \theta, obtained by preprocessing the
% t-th subresultant matrix.
[t, lambda, mu, rho, theta, gm_fx, gm_gx, gm_hx] = ...
    Get_GCD_Degree_3Polys(fx, gx, hx, limits_t, rank_range);






% If the computed degree of the GCD is 0 then the polynomials f(x), g(x)
% and h(x) are coprime.
if t == 0 
    
    fprintf([mfilename ' : ' ...
        sprintf('Polynomials f(x), g(x) and h(x) are coprime \n')])
    
    % Set polynomials d(x), u(x), v(x) and w(x)
    dx_o = 1;
    ux_o = fx;
    vx_o = gx;
    wx_o = hx;
    
    % Set optimal values lambda, rho and mu
    lambda_o = 1;
    mu_o = 1; 
    rho_o = 1;
    
    return
    
end







% If finding the GCD fails, set the degree of the GCD to be 1.
if isempty(t)
    t = 1;
end


% Normalise f(x) and g(x) by geometric mean to obtain fx_n and gx_n.
% Normalise by geometric mean obtained by entries of f(x) and g(x) in the
% subresultant S_{t}
fx_n = fx ./ gm_fx;
gx_n = gx ./ gm_gx;
hx_n = hx ./ gm_hx;



% % Get the optimal column of the sylvester matrix to be removed. Where
% % removal of the optimal column gives the minmal residual in (Ak x = ck)



% Get the set of preprocessed polynomials  
%   \lambda_{t}  \tilde{f}_{t}(\omega) 
%   \mu_{t} \tilde{g}_{t}(\omega)
%   \rho_{t} \tilde{h}_{t}(\omega)
lambda_fw = lambda .* GetWithThetas(fx_n, theta);
mu_gw  = mu .* GetWithThetas(gx_n, theta);
rho_hw = rho .* GetWithThetas(hx_n, theta);


% Build the t-th preprocessed subresultant matrix S_{t}(f, g, h)
St_preproc = BuildSubresultant_3Polys(lambda_fw, mu_gw, rho_hw, t);




% Get index of optimal column for removal from the t-th preprocessed
% subresultant matrix
[~, idx_col] = GetMinDistance(St_preproc);





% % Get low rank approximation of the t-th subresultant matrix S_{t}(f, g,
% h) and the coefficients of the cofactor polynomials u(x), v(x) and w(x)
[fx_lr, gx_lr, hx_lr, ux_lr, vx_lr, wx_lr, lambda_lr, mu_lr, ...
    rho_lr, theta_lr] = LowRankApproximation_3Polys(fx_n, ...
    gx_n, hx_n, lambda, mu, rho, theta, t, idx_col);




% Get the coefficients of the GCD d(x) by APF or other method.
[fx_alr, gx_alr, hx_alr, dx_alr, ux_alr, vx_alr, wx_alr, lambda_alr, ...
    mu_alr, rho_alr, theta_alr] = ...
    APF_3Polys(fx_lr, gx_lr, hx_lr, ux_lr, vx_lr, wx_lr, ...
    lambda_lr, mu_lr, rho_lr, theta_lr, t);






% % Set function outputs

% Output vectors of coefficients of the polynomials f(x), g(x) and h(x)
fx_o = fx_alr;
gx_o = gx_alr;
hx_o = hx_alr;

% Output the coefficients of the GCD
dx_o = dx_alr;

% Output the coefficients of the cofactor polynomials u(x), v(x) and w(x)
ux_o = ux_alr;
vx_o = vx_alr;
wx_o = wx_alr;

% Output the Optimal values \lambda_{t}, \mu_{t}, \rho_{t} and \theta_{t}
lambda_o = lambda_alr;
mu_o = mu_alr;
rho_o = rho_alr;
theta_o = theta_alr;




end









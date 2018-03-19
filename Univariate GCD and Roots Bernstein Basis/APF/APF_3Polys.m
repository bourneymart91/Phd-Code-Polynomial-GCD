function [fx_lr, gx_lr, hx_lr, dx_lr, ux_lr, vx_lr, wx_lr, ...
    lambda_lr, mu_lr, rho_lr, theta_lr] = ...
    APF_3Polys(fx, gx, hx, ux, vx, wx, lambda, mu, rho, theta, k)
%
% % Inputs
%
% fx : (Vector) The coefficients of the polynomial f(x)
%
% gx : (Vector) The coefficients of the polynomial g(x)
%
% hx : (Vector) The coefficients of the polynomial h(x)
%
% ux : (Vector) The coefficients of the polynomial u(x)
%
% vx : (Vector) The coefficients of the polynomial v(x)
%
% wx : (Vector) The coefficients of the polynomial w(x)
%
% lambda : (Float) Optimal value of \lambda_{k}
%
% mu : (Float) Optimal value of \mu_{k}
%
% rho : (Float) Optimal value of \rho_{k}
%
% theta : (Float) Optimal value of \theta_{k}
% 
% k : (Int) The degree of polynomial d_{k}(x)
%
% % Outputs
%
% fx_lr : (Vector) The coefficients of the output polynomial f(x)
%
% gx_lr : (Vector) The coefficients of the output polynomial g(x)
%
% hx_lr : (Vector) The coefficients of the output polynomial h(x)
%
% dx_lr : (Vector) The coefficients of the output polynomial d(x)
%
% ux_lr : (Vector) The coefficients of the output polynomial u(x)
%
% vx_lr : (Vector) The coefficients of the output polynomial v(x)
%
% wx_lr : (Vector) The coefficients of the output polynomial w(x)
%
% lambda_lr : (Float) Optimal value of \lambda_{k}
%
% mu_lr : (Float) Optimal value of \mu_{k}
% 
% rho_lr : (Float) Optimal value of \rho_{k}
%
% theta_lr : (Float) Optimal value of \theta_{k}





% if not 11 input arguments
if (nargin ~= 11)
    error('Not enough input arguments');
end



global SETTINGS;


switch SETTINGS.APF_METHOD

    case 'Standard APF Nonlinear'
        
        error([mfilename ' : Code Not Yet Developed']);
        
        % This function does not yet exist 
        
        [fx_lr, gx_lr, hx_lr, dx_lr, ux_lr, vx_lr, wx_lr, lambda_lr, mu_lr, rho_lr, theta_lr] = ...
            APF_Nonlinear_3Polys(fx, gx, hx, ux, vx, wx, lambda, mu, rho, theta, k);
 

    case 'Standard APF Linear'

        error([mfilename ' : Code Not Yet Developed']);
        
        % Get f(\omega), g(\omega) and h(\omega)
        lambda_fw = lambda .* GetWithThetas(fx1, theta);
        mu_gw = mu .* GetWithThetas(gx, theta);
        rho_hw = GetWithThetas(hx, theta);
        
        % Get u(\omega), v(\omega) and w(\omega)
        uw = GetWithThetas(ux, theta);
        vw = GetWithThetas(vx, theta);
        ww = GetWithThetas(wx, theta);
        
        % Get APF and d(\omega)
        % This function does not yet exist
        [lambda_fw_lr, mu_gw_lr, rho_dw_lr, uw_lr, vw_lr, ww_lr] = ...
            APF_Linear_3Polys(lambda_fw, mu_gw, rho_hw, uw, vw, ww, k);
        
        % Get f(x) and g(x) after linear APF function
        fx_lr = GetWithoutThetas(lambda_fw_lr, theta) ./ lambda;
        gx_lr = GetWithoutThetas(mu_gw_lr, theta) ./ mu;
        hx_lr = GetWithoutThetas(rho_hw_lr, theta) ./ rho;
        
        % Get u(x) and v(x) after linear APF function
        ux_lr = GetWithoutThetas(uw_lr, theta);
        vx_lr = GetWithoutThetas(vw_lr, theta);
        wx_lr = GetWithoutThetas(ww_lr, theta);
        
        % Get d(x) after linear APF function
        dx_lr = GetWithoutThetas(dw_lr, theta);
        
        
        lambda_lr = lambda;
        mu_lr = mu;
        rho_lr = rho;
        theta_lr = theta;
        
    case 'None'
        
        fx_lr = fx;
        gx_lr = gx;
        hx_lr = hx;
        
        ux_lr = ux;
        vx_lr = vx;
        wx_lr = wx;
        
        % Get the 
        lambda_fw   = lambda    .* GetWithThetas(fx, theta);
        mu_gw       = mu        .* GetWithThetas(gx, theta);
        rho_hw      = rho       .* GetWithThetas(hx, theta);
        
        uw = GetWithThetas(ux, theta);
        vw = GetWithThetas(vx, theta);
        ww = GetWithThetas(wx, theta);
        
        % Get the vector of coefficients of the polynomial d(\omega)
        dw_lr = GetGCDCoefficients_3Polys(uw, vw, ww, ...
            lambda_fw, mu_gw, rho_hw, k);
        
        % Get the vector of coefficients of the polynomial d(x)
        dx_lr = GetWithoutThetas(dw_lr, theta);
        
        
        % lambda, mu, rho and theta are unaltered from input
        lambda_lr   = lambda;
        mu_lr       = mu;
        rho_lr      = rho;
        theta_lr    = theta;
        
        % No iterations required in this method
        SETTINGS.APF_REQ_ITE = 0;
        
    otherwise
        error([mfilename sprintf(' : Error : %s is not a valide APF Method',SETTINGS.APF_METHOD)])
end


end
function [ux_lra, vx_lra, fx_lra, gx_lra, dx_lra, alpha_lra, theta_lra] ...
    = APF_2Polys(ux, vx, fx, gx, alpha, theta, t)
% Get the coefficients of the GCD d(x), of f(x) and g(x), given u(x) and v(x).
% This function has two branches. d(x) can be computed either by utilising
% [u(x), v(x),f(x) and g(x)], or just [u(x) and f(x)]
%
% % Inputs.
%
% ux : (Vector) Coefficients of polynomial u(x)
%
% vx : (Vector) Coefficients of polynomial v(x)
%
% fx : (Vector) Coefficients of polynomial f(x) 
%
% gx : (Vector) Coefficients of polynomial g(x)
%
% alpha : (Float) Optimal value of \alpha
%
% theta : (Float) Optimal value of \theta
%
% k : (Int) Degree of the GCD d(x)

global SETTINGS
switch SETTINGS.APF_METHOD
    case 'None'
        
        % Get f(\omega) and g(\omega)
        fw = GetWithThetas(fx,theta);
        a_gw = alpha.* GetWithThetas(gx,theta);
        
        % Get u(\omega) and v(\omega)
        uw = GetWithThetas(ux,theta);
        vw = GetWithThetas(vx,theta);
        
        % Get d(\omega)
        [dw] = GetGCDCoefficients_2Polys(uw, vw, fw, a_gw,t);
        
        % Get d(x) output
        dx_lra = GetWithoutThetas(dw,theta);
        
        % Get u(x) and v(x) output unchanged from input
        ux_lra = ux;
        vx_lra = vx;
        
        % Get f(x) and g(x) output, unchanged from input
        fx_lra = fx;
        gx_lra = gx;
        
        % Get \alpha and \theta output, unchanged from input
        alpha_lra = alpha;
        theta_lra = theta;
        
        SETTINGS.APF_REQ_ITE = 0;
        
    case 'Standard APF Nonlinear'
        
        error([mfilename ' : ' 'Code not completed']);
        SETTINGS.APF_REQ_ITE = 0;
        
    case 'Standard APF Linear'
        
        error([mfilename ' : ' 'Code not completed']);
        SETTINGS.APF_REQ_ITE = 0;
        
    case 'Last Non-zero Row'
        
        m = GetDegree(fx);
        n = GetDegree(gx);
        
        
        fw = GetWithThetas(fx, theta);
        a_gw = alpha.* GetWithThetas(gx, theta);
     
        S1 = BuildT(fw, a_gw, 1)';
        
        [~,r] = qr(S1);
            
        columns = m + n - t  : m + n;
        
        dw = r(m + n - t, columns)';
        dw = dw./dw(1);
        dx_lra = GetWithoutThetas(dw, theta);
        
        
        fx_lra = fx;
        gx_lra = gx;
        alpha_lra = alpha;
        theta_lra = theta;
        
        ux_lra = ux;
        vx_lra = vx;
        
    otherwise
        error('err')
end



end

function [fx_lr, gx_lr, dx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = ...
    APF_2Polys(fx, gx, ux, vx, alpha, theta, t)
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x) 
%
% gx : (Vector) Coefficients of polynomial g(x) 
%
% ux : (Vector) Coefficients of polynomial u(x) 
%
% vx : (Vector) Coefficients of the polynomial v(x)
% 
% alpha : (Float) Optimal value of \alpha
%
% theta : (Float) Optimal value of \theta
% 
% t : (Int) Degree of polynomial d(x)
%
% % Outputs
%
% fx_lr : (Vector) Coefficients of the polynomial f(x)
%
% gx_lr : (Vector) Coefficients of the polynomial g(x)
%
% dx_lr : (vector) Coefficients of the polynomial d(x)
%
% ux_lr : (vector) Coefficients of the polynomial u(x)
%
% vx_lr : (Vector) Coefficients of the polynomial v(x)





global SETTINGS;


switch SETTINGS.APF_METHOD

    % Standard APF Nonlinear
    % Standard APF Linear
    % None
    
    case 'Standard APF Nonlinear'
        
        [fx_lr, gx_lr, dx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = ...
            APF_Nonlinear_2Polys(fx, gx, ux, vx, alpha, theta, t);
 
        res1 = GetResidual(fx, gx, ux, vx, dx);
        
        fprintf('Residual : %e', res1)
        

    case 'Standard APF Linear'

        % Get f(\omega) and g(\omega)
        fw = GetWithThetas(fx, theta);
        a_gw = alpha.*GetWithThetas(gx, theta);
        
        % Get u(\omega) and v(\omega)
        uw = GetWithThetas(ux, theta);
        vw = GetWithThetas(vx, theta);
        
        % Get APF and d(\omega)
        [fw_lr, a_gw_lr, dw_lr, uw_lr, vw_lr] = APF_Linear_2Polys(fw, a_gw, uw, vw, t);
        
        % Get f(x) and g(x) after linear APF function
        fx_lr = GetWithoutThetas(fw_lr, theta);
        gx_lr = GetWithoutThetas(a_gw_lr, theta) ./ alpha;
        
        % Get u(x) and v(x) after linear APF function
        ux_lr = GetWithoutThetas(uw_lr, theta);
        vx_lr = GetWithoutThetas(vw_lr, theta);
        
        % Get d(x) after linear APF function
        dx_lr = GetWithoutThetas(dw_lr, theta);
        
        alpha_lr = alpha;
        theta_lr = theta;
     
        
    case 'Last Non-zero Row'
        
        m = GetDegree(fx);
        n = GetDegree(gx);
        
        
        fw = GetWithThetas(fx, theta);
        a_gw = alpha.* GetWithThetas(gx, theta);
     
        S1 = BuildSubresultant_2Polys(fw, a_gw, 1)';
        
        [~,r] = qr(S1);
            
        columns = m + n - t  : m + n;
        
        dx_lr = r(m + n - t, columns)';
        
        
        %dx_lr = GetWithoutBinomials(dx_lr);
        %dx_lr = GetWithoutThetas(dx_lr, theta);
        
        dx_lr = dx_lr./dx_lr(1);

        
        
        fx_lr = fx;
        gx_lr = gx;
        alpha_lr = alpha;
        theta_lr = theta;
        
        ux_lr = ux;
        vx_lr = vx;    
        
        
    case 'None'
        
     
        
        uw = GetWithThetas(ux, theta);
        vw = GetWithThetas(vx, theta);
        fw = GetWithThetas(fx, theta);
        a_gw = alpha .* GetWithThetas(gx, theta);
        
        dw_lr = GetGCDCoefficients_2Polys(uw, vw, fw, a_gw, t);
        
        alpha_lr = alpha;
        theta_lr = theta;
        
        fx_lr = fx;
        gx_lr = gx;
        ux_lr = ux;
        vx_lr = vx;
        dx_lr = GetWithoutThetas(dw_lr, theta);
        
        a_gx = alpha.* gx;
        
        res1 = GetResidual(fx, a_gx, ux, vx, dx_lr);
        res2 = GetResidual(fw, a_gw, uw, vw, dw_lr);
        
        LineBreakLarge()
        fprintf('Computing GCD Coefficients by APF \n')
        fprintf('Residual f(x), a_g(x), u(x), v(x) : %e \n', res1)
        fprintf('Residual f(w), a_g(w), u(w), v(w) : %e \n', res2)
        
        
        SETTINGS.APF_REQ_ITE = 0;
        
    otherwise
        error([mfilename sprintf(' : Error : %s is not a valid APF Method',SETTINGS.APF_METHOD)])
end


end


function residual = GetResidual(fx, gx, ux, vx, dx)

m = GetDegree(fx);
n = GetDegree(gx);
t = GetDegree(dx);

HCG = BuildHCG_2Polys(ux,vx,t);
vec_fg = [fx;gx];

residual = norm(vec_fg - HCG*dx);






end

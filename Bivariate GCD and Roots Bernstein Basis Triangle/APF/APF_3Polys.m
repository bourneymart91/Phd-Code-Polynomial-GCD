function [fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, dxy_lr, ...
    lambda_lr, mu_lr, rho_lr, th1_lr, th2_lr] = APF_3Polys(fxy, gxy, hxy, uxy, ...
    vxy, wxy, lambda, mu, rho, th1, th2, m, n, o, k)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of f(x,y)
%
% gxy : (Matrix) Coefficients of g(x,y)
%
% hxy : (Matrix) Coefficients of h(x,y)
%
% uxy : (Matrix) Coefficients of u(x,y)
%
% vxy : (Matrix) Coefficients of v(x,y)
%
% wxy : (Matrix) Coefficietns of w(x,y)
%
% alpha : (Float) Optimal value of alpha
%
% beta : (Float) Optimal value of beta
%
% gamma : (Float)
%
% th1 : (Float) Optimal value of theta_{1}
%
% th2 : (Float) Optimal value of theta_{2}
% 
% m : (Int) Total degree of polynomial f(x,y)
%
% n : (Int) Total degree of polynomial g(x,y)
%
% o : (Int) Total degree of polynomial h(x,y)
%
% k : (Int)
%
%
% % Outputs
%
% [fxy_lr, gxy_lr, hxy_lr] : (Matrix)(Matrix)(Matrix)
% 
% [uxy_lr, vxy_lr, wxy_lr] : (Matrix)(Matrix)(Matrix)
%
% dxy_lr : (Matrix)
% 
% alpha_lr : (Float)
%
% [th1_lr, th2_lr] : (Float)(Float)

global SETTINGS

switch SETTINGS.APF_METHOD
    case 'Standard APF Nonlinear'
        
        error([mfilename ' : This APF method is not yet developed.'])
        SETTINGS.LOW_RANK_APPROX_REQ_ITE = 0;
        
    case 'Standard APF Linear'
        
        error([mfilename ' : This APF method is not yet developed. '])
        SETTINGS.LOW_RANK_APPROX_REQ_ITE = 0;
        
    case 'None'
        
        % Get f(\omega_{1},\omega_{2}) and g(\omega_{1},\omega_{2})
        fww = GetWithThetas(fxy, m, th1, th2);        
        gww = GetWithThetas(gxy, n, th1, th2);
        hww = GetWithThetas(hxy, o, th1, th2);
        
        % Get u(\omega_{1},\omega_{2}) and v(\omega_{1},\omega_{2})
        uww = GetWithThetas(uxy, m - k, th1, th2);
        vww = GetWithThetas(vxy, n - k, th1, th2);
        www = GetWithThetas(wxy, o - k, th1, th2);
        
        % Get d(\omega_{1},\omega_{2})
        dww = GetGCDCoefficients_3Polys(...
            lambda .* fww, ...
            mu .* gww, ...
            rho .* hww, ...
            uww, vww, www, ...
            m, n, o, k);
        
        % Get f(x,y) and g(x,y)
        fxy_lr = fxy;
        gxy_lr = gxy;
        hxy_lr = hxy;
        
        % Get u(x,y) and v(x,y)
        uxy_lr = uxy;
        vxy_lr = vxy;
        wxy_lr = wxy;
        
        % Get \alpha, \theta_{1} and \theta_{2}
        lambda_lr = lambda;
        mu_lr = mu;
        rho_lr = rho;
        
        th1_lr = th1;
        th2_lr = th2;
        
        dxy_lr = GetWithoutThetas(dww, k, th1, th2);
        
        SETTINGS.APF_REQ_ITE = 0;
        
    otherwise
        error([mfilename ' : ' SETTINGS.APF_METHOD ': Not a valid APF method.'])
end

end
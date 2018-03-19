function [fxy_lr, gxy_lr, uxy_lr, vxy_lr, dxy_lr, alpha_lr, th1_lr, th2_lr] = ...
    APF_2Polys(fxy, gxy, uxy, vxy, alpha, th1, th2, m, n, k)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% uxy : (Matrix) Coefficients of the polynomial u(x,y), the quotient
% polynomial such that f(x,y)/u(x,y) = d(x,y)
%
% vxy : (Matrix) Coefficients of the polynomial u(x,y), the quotient
% polynomial such that g(x,y)/v(x,y) = d(x,y)
%
% alpha : (Float)
%
% th1 : (Float)
%
% th2 : (Float)
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% k : Total degree of d_{k}(x,y)
%
%
% % Outputs
%
% fxy_lr : 
%
% gxy_lr :
%
% uxy_lr : (Matrix) Coefficients of the polynomial u(x,y), the quotient
% polynomial such that f(x,y)/u(x,y) = d(x,y)
%
% vxy_lr : (Matrix) Coefficients of the polynomial u(x,y), the quotient
% polynomial such that g(x,y)/v(x,y) = d(x,y)
%
% dxy_lr : (Matrix) Coefficients of the polynomial d(x,y)
%
% alpha_lr : (Float) 
%
% th1_lr : (Float)
%
% th2_lr : (Float)


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
        a_gww = alpha .* GetWithThetas(gxy, n, th1, th2);
        
        % Get u(\omega_{1},\omega_{2}) and v(\omega_{1},\omega_{2})
        uww = GetWithThetas(uxy, m-k, th1, th2);
        vww = GetWithThetas(vxy, n-k, th1, th2);
   
        % Get d(\omega_{1},\omega_{2})
        dww = GetGCDCoefficients_2Polys(fww,a_gww,uww,vww,m,n,k);
        
        % Get f(x,y) and g(x,y)
        fxy_lr = fxy;
        gxy_lr = gxy;
        
        % Get u(x,y) and v(x,y)
        uxy_lr = uxy;
        vxy_lr = vxy;
        
        % Get \alpha, \theta_{1} and \theta_{2}
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        dxy_lr = GetWithoutThetas(dww,k,th1,th2);
        
        SETTINGS.APF_REQ_ITE = 0;
        
    otherwise
        error([mfilename ' : ' SETTINGS.APF_METHOD ': Not a valid APF method.'])
end

end
function [fxy_lr, gxy_lr, uxy_lr, vxy_lr, dxy_lr, alpha_lr, th1_lr, th2_lr ] ...
    = APF_Bivariate_2Polys(fxy, gxy, uxy, vxy, m, n, t, alpha, th1, th2)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% uxy : (Matrix) Coefficients of cofactor polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of cofactor polynomial v(x,y)
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% t : (Int) Total degree of d(x,y)
%
% alpha : (Float) : Opitmal value of \alpha
%
% th1 : (Float) : Optimal value of \theta_{1}
%
% th2 : (Float) : Optimal value of \theta_{2}

global SETTINGS



switch SETTINGS.APF_METHOD
    case 'None'
        
        % Get f(w,w) from f(x,y) and g(w,w) from g(x,y)
        fww = GetWithThetas(fxy, th1, th2);
        gww = GetWithThetas(gxy, th1, th2);
        
        % Get u(w1,w2) and v(w1,w2) from u(x,y) and v(x,y)
        uww = GetWithThetas(uxy, th1, th2);
        vww = GetWithThetas(vxy, th1, th2);
        
        % Get the coefficients of the GCD d(x,y)
        [dww_lr] = GetGCDCoefficients_Bivariate_2Polys(fww, alpha.*gww, uww, vww, m, n, t);
        
        % Get d(x,y) from d(w1,w2)
        dxy_lr = GetWithoutThetas(dww_lr, th1, th2);
        
        % Get outputs
        fxy_lr = fxy;
        gxy_lr = gxy;
        
        uxy_lr = uxy;
        vxy_lr = vxy;
        
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
    case 'Standard Nonlinear APF'
        
        error([mfilename ' : Method not developed']);
        
    case 'Standard Linear APF'
        
        error([mfilename ' : Method not developed']);
        
    otherwise
        
        error([mfilename ' : Error APF Method is either None, Standard Nonlinear APF or Standard Linear APF '])
        
end


end
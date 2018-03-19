function [fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, dxy_lr, alpha_lr, th1_lr, th2_lr ] ...
    = APF_Bivariate_3Polys(fxy, gxy, hxy, uxy, vxy, wxy, m, n, o, t, alpha, th1, th2)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of polynomial f(x,y)
%
% uxy : (Matrix) Coefficients of cofactor polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of cofactor polynomial v(x,y)
% 
% wxy : (Matrix) Coefficients of cofactor polynomial w(x,y)
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% o : (Int) Total degree of h(x,y)
%
% t : (Int) Total degree of d(x,y)
%
% alpha : (Float) : Optimal value of \alpha
%
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}
%
% % Outputs
%
% fxy_lr : (Matrix)
%
% gxy_lr : (Matrix) 
%
% hxy_lr : (Matrix) 
%
% uxy_lr : (Matrix)
%
% vxy_lr : (Matrix) 
%
% wxy_lr : (Matrix)
%
% dxy_lr : (Matrix)
%
% alpha_lr : (Float)
%
% th1_lr : (Float)
%
% th2_lr : (Float)


global SETTINGS



switch SETTINGS.APF_METHOD
    case 'None'
        
        % Get preprocessed form of f(x,y), g(x,y) and h(x,y)
        fww = GetWithThetas(fxy, th1, th2);
        gww = GetWithThetas(gxy, th1, th2);
        hww = GetWithThetas(hxy, th1, th2);
        
        % Get preprocessed form of u(x,y), v(x,y) and w(x,y)
        uww = GetWithThetas(uxy, th1, th2);
        vww = GetWithThetas(vxy, th1, th2);
        www = GetWithThetas(wxy, th1, th2);
        
        % Get the coefficients of the GCD d(x,y)
        [dww_lr] = GetGCDCoefficients_Bivariate_3Polys(fww, alpha.*gww, hww, uww, vww, www, m, n, o, t);
        
        % Get unprocessed form d(x,y)
        dxy_lr = GetWithoutThetas(dww_lr, th1, th2);
        
        % Set function outputs
        fxy_lr = fxy;
        gxy_lr = gxy;
        hxy_lr = hxy;
        
        uxy_lr = uxy;
        vxy_lr = vxy;
        wxy_lr = wxy;
        
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
    case 'Standard Nonlinear APF'
        
        error([mfilename ' : Method not developed']);
        
    case 'Standard Linear APF'
        
        error([mfilename ' : Method not developed']);
        
    otherwise 
        
        error([mfilename ' : Error'])
        
end


end
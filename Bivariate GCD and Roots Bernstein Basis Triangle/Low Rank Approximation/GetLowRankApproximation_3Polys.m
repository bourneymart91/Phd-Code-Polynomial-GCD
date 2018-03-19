function [fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr,...
    alpha_lr, beta_lr, gamma_lr, th1_lr, th2_lr] = ...
    GetLowRankApproximation_3Polys(fxy, gxy, hxy, alpha, beta, ...
    gamma, th1, th2, m, n, o, k)
% Get the low rank approximation of the Sylvester matrix S_{t}(f,g), and
% return the polynomials f_lr(x,y) and g_lr(x,y) which are the perturbed forms of
% input polynomials f(x,y) and g(x,y).
% lr = low rank.
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
% 
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of polynomial h(x,y)
%
% alpha : (Float) Optimal value of \alpha
%
% beta : (Float) Optimal value of \beta
%
% th1 : (Float) \theta_{1}
%
% th2 : (Float) \theta_{2}
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% o : (Int) Total degree of h(x,y)
%
% % Outputs
%
% fxy_lr : (Matrix) 
%
% gxy_lr : (Matrix) 
%
% hxy_lr : (Matrix)  Coefficients of polynomial f_lr(x,y) which is 
% used in the low rank approximation of S_{t}(f,g).
%
% uxy_lr : (Matrix) 
%
% vxy_lr : (Matrix) 
%
% wxy_lr : (Matrix) 
%
% alpha : (Float) \alpha
%
% beta : (Float)
%
% gamma : (Float)
%
% th1 : (Float) 
%
% th2 : (Float) \theta_{2}

global SETTINGS

switch SETTINGS.LOW_RANK_APPROX_METHOD
    
    
    case 'None'
        
        % Get f(\omega_{1},\omega_{2}), g(\omega_{1},\omega_{2}) and h(\omega_{1},\omega_{2})
        fww = GetWithThetas(fxy, m, th1, th2);
        gww = GetWithThetas(gxy, n, th1, th2);
        hww = GetWithThetas(hxy, o, th1, th2);
               
        % Get polynomials u(\omega_{1},\omega_{2}) and
        % v(\omega_{1},\omega_{2}).
        [uww, vww, www] = GetCofactors_3Polys(...
            alpha .* fww,...
            beta .* gww,...
            gamma .* hww, ...
            m, n, o, k);
        
        fxy_lr = fxy;
        gxy_lr = gxy;
        hxy_lr = hxy;
                
        % Get u(x,y), v(x,y) and w(x,y)
        uxy_lr = GetWithoutThetas(uww, m - k, th1, th2);
        vxy_lr = GetWithoutThetas(vww, n - k, th1, th2);
        wxy_lr = GetWithoutThetas(www, o - k, th1, th2);
        
        alpha_lr = alpha;
        beta_lr = beta;
        gamma_lr = gamma;
        
        th1_lr = th1;
        th2_lr = th2;
        
        SETTINGS.LOW_RANK_APPROX_REQ_ITE = 0;
        
        
    case 'Standard STLN'
        
        error('Not Yet Complete')
        
       
        
    case 'Standard SNTLN'
        
        error('Not Yet Complete')
       
end
end
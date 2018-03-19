function [fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr,th2_lr] = ...
    LowRankApprox_Bivariate_2Polys(fxy, gxy, alpha, th1, th2, m, n, k, k1, k2, idx_col)
% Compute low rank approximation of the sylvester matrix formed from
% coefficients of f(x,y) and g(x,y). Return the modified forms of f(x,y)
% and g(x,y).
%
% % Inputs
% 
% [fxy, gxy] : Coefficients of polynomials f(x,y) and g(x,y)
%
% alpha : \alpha
%
% th1 : \theta_{1}
%
% th2 : \theta_{2}
%
% [m, n] : Total degree of f(x,y) and g(x,y)
%
% k : Total degree of d(x,y)
% 
% k1 : Degree of d(x,y) with respect to x
%
% k2 : Degree of d(x,y) with respect to y
%
% idx_col : Index of optimal column for removal from S_{}(f,g)
%
% % Outputs.
%
% fxy_lr : Coefficients of polynomial f(x,y) with added perturbations
% 
% gxy_lr : Coefficients of polynomial g(x,y) with added perturbations
%
% alpha : Refined alpha
% 
% th1 : Refined theta_{1}
% 
% th2 : Refined theta_{2}


% Global Settings
global SETTINGS

% Note. Use the notation 'lr' for low rank approximation.

switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    
    case 'Standard STLN'
        
        % Get preprocessed polynomials
        fww = GetWithThetas(fxy, th1, th2);
        gww = GetWithThetas(gxy, th1, th2);
        
        % Get Low rank approximation by STLN
        [fww_lr, a_gww_lr, uww_lr, vww_lr] = STLN(fww, alpha.*gww, m, n, k, k1, k2, idx_col);
        
        % Remove alpha from alpha.*g(w,w)
        gww_lr = a_gww_lr ./ alpha;
        
        % Get f(x,y) from f(w,w).
        % Get g(x,y)_lr from g(w,w).
        fxy_lr = GetWithoutThetas(fww_lr, th1, th2);
        gxy_lr = GetWithoutThetas(gww_lr, th1, th2);
        
        % Get u(x,y) from u(w,w)
        % Get u(x,y) from u(w,w)
        uxy_lr = GetWithoutThetas(uww_lr, th1, th2);
        vxy_lr = GetWithoutThetas(vww_lr, th1, th2);
        
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        % Get distance between f(x,y) and f(x,y) low rank approximation by
        % STLN.
        test = (fxy_lr - fxy);
        norm(test)
        
        % Get distance between f(w,w) and f(w,w) low rank approximation by
        % STLN.
        test2 = (fww_lr - fww);
        norm(test2)
                
        % Note that minimial perturbations of f(w,w) found by low rank approximation
        % may be small, but CAN be large perturbations of f(x,y).
        
        
        
        
    case 'Standard SNTLN'
        
        % Get low rank approximation by SNTLN
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
                    SNTLN(fxy, gxy, alpha, th1, th2, m, n, k, k1, k2, idx_col);
        
        
        % Update f(x,y) g(x,y) theta1 theta2 and alpha to their new values post
        % SNTLN
        
        fww = GetWithThetas(fxy, th1, th2);
        gww = GetWithThetas(gxy, th1, th2);
        
        % Get f(w,w)_lr from f(x,y)_lr
        fww_lr = GetWithThetas(fxy, th1_lr, th2_lr);
        
        % Get g(w,w)_lr from g(x,y)_lr
        gww_lr = GetWithThetas(gxy, th1_lr, th2_lr);
        
        % 
        test1 = (fww - fww_lr);
        norm(test1)
        test2 = (fxy - fxy_lr);
        norm(test2)

        
    case 'None'
        
        fww = GetWithThetas(fxy, th1, th2);
        gww = GetWithThetas(gxy, th1, th2);
        a_gww = alpha.*gww; 
        
        [uww_lr,vww_lr] = GetQuotients_Bivariate_2Polys(fww, a_gww, m, n, k, k1, k2);
        
        uxy_lr = GetWithoutThetas(uww_lr, th1, th2);
        vxy_lr = GetWithoutThetas(vww_lr, th1, th2);
        
        fxy_lr = fxy;
        gxy_lr = gxy;
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        fww_lr = fww;
        gww_lr = gww;
        
        SETTINGS.LOW_RANK_APPROX_REQ_ITE = 0;
        
    otherwise
        error([mfilename ' : ' 'LOW_RANK_APPROXIMATION_METHOD' ...
            'must be set to either (Standard STLN) or (Standard SNTLN) or (None)'])
end

        % Build Sylvester matrix of input polynomials.
        S1 = BuildT_Bivariate_2Polys(fxy, gxy, m, n, k, k1, k2);
        
        % Build Sylvester matrix of low rank approx of input polys.
        S2 = BuildT_Bivariate_2Polys(fxy_lr, gxy_lr, m, n, k, k1, k2);
        
        % Build Sylvester matrix of polys in preprocessed form
        S3 = BuildT_Bivariate_2Polys(fww, alpha.*gww, m, n, k, k1, k2);
        
        % Build Sylvester matrix of polys in preprocessed form with pert.
        S4 = BuildT_Bivariate_2Polys(fww_lr, alpha.*gww_lr, m, n, k, k1, k2);
        
        % Get Singular Values
        vec_SingularValues_1 = svd(S1);
        vec_SingularValues_2 = svd(S2);
        vec_SingularValues_3 = svd(S3);
        vec_SingularValues_4 = svd(S4);
        
        
        % Plot singular values
        figure_name = sprintf([mfilename ' : ' 'Singular Values']);
        figure('name', figure_name )
        plot(log10(vec_SingularValues_1),'-s', 'DisplayName','f(x,y) g(x,y)')
        hold on
        plot(log10(vec_SingularValues_2),'-*', 'DisplayName','f(x,y) g(x,y) low rank')
        plot(log10(vec_SingularValues_3),'-s', 'DisplayName','f(w,w) g(w,w)')
        plot(log10(vec_SingularValues_4),'-s', 'DisplayName','f(w,w) g(w,w) low rank')
        legend(gca,'show');
        hold off

end
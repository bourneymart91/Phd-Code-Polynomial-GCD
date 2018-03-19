function [fx_lr, gx_lr, hx_lr, ux_lr, vx_lr, wx_lr, lambda_lr,...
    mu_lr, rho_lr, theta_lr] = ...
    LowRankApproximation_3Polys(fx, gx, hx, lambda, mu, rho, theta, k, idx_col)
% Get the low rank approximation of the Sylvester subresultant matrix
% S_{k}(f,g)
%
% % Inputs.
%
% fx : (Vector) Vector containing the coefficients of the polynomial f(x)
%
% gx : (Vector) Vector containing the coefficients of the polynomial g(x)
%
% hx : (Vector) Vector containing the coefficients of the polynomial h(x)
%
% alpha : (Float) Optimal value of \lambda_{k}
%
% beta : (Float) Optimal value of \mu_{k}
%
% theta : (Float) Optimal value of \theta_{k}
%
% k : (Int) Degree of GCD of f(x) and g(x)
%
% idx_col : (Int) Index of optimal column for removal from S_{k}(f,g)
%
% % Outputs
%
% fx_lr : (Vector) Vector containing the coefficients of the polynomial f(x)
%
% gx_lr : (Vector) Vector containing the coefficients of the polynomial g(x)
%
% hx_lr : (Vector) Vector containing the coefficients of the polynomial h(x)
%
% ux_lr : (Vector) Vector containing the coefficients of the polynomial u(x)
%
% vx_lr : (Vector) Vector containing the coefficients of the polynomial v(x)
%
% wx_lr : (Vector) Vector containing the coefficients of the polynomial w(x)
%
% lambda_lr : (Float) optimal value of \lambda_{k}
%
% mu_lr : (Float) Optimal value of \mu_{k}
%
% rho_lr : (Float) Optimal value of \rho_{k}
%
% theta_lr : (Float) Optimal value of \theta_{k}

global SETTINGS


% Nine input arguments expected
if (nargin ~= 9)
   error('Not enough input arguments') 
end


switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    
    case 'Standard STLN'
        
        error([mfilename ' : Code Not Yet Completed'])
        
        % The file STLN_3Polys does not yet exist!!!
        
        % Get f(\omega) and g(\omega)
        lambda_fw   = lambda 	.* GetWithThetas(fx, theta);
        mu_gw       = mu        .* GetWithThetas(gx, theta);
        rho_hw      = rho       .* GetWithThetas(hx, theta);
        
        % Performe STLN to get low rank approximation of S(f,g)
        [lambda_fw_LR, mu_gw_LR, rho_hw_LR, uw_LR, vw_LR, ww_LR] = ...
            STLN_3Polys(lambda_fw, mu_gw, rho_hw, k, idx_col);
        
        
        
        % Get f(x), g(x) and h(x) from low rank approximation.
        fx_lr = GetWithoutThetas(lambda_fw_LR, theta)   ./ lambda;
        gx_lr = GetWithoutThetas(mu_gw_LR, theta)       ./ mu;
        hx_lr = GetWithoutThetas(rho_hw_LR, theta)      ./ rho;
        
        % Get u(x), v(x) and w(x) from low rank approximation
        ux_lr = GetWithoutThetas(uw_LR, theta);
        vx_lr = GetWithoutThetas(vw_LR, theta);
        wx_lr = GetWithoutThetas(ww_LR, theta);
        
        % \lamda_{k}, \rho_{k}, \mu_{k} and \theta_{k} are unchanged by 
        % STLN
        lambda_lr   = lambda;
        mu_lr       = mu;
        rho_lr      = rho;
        theta_lr = theta;

        % Graph Plotting
        Plot_LowRank_SingularValues(fx, gx, hx, fx_lr, gx_lr, hx_lr, ...
            lambda_fw, mu_gw, rho_hw, lambda_fw_LR, mu_gw_LR, rho_hw_LR, k)

        
        
        
        
        
        
    case 'Standard SNTLN' % Structured Non-Linear Total Least Norm
        
        error([mfilename ' : Code Not Yet Completed'])
        
        % Perform SNTLN to get low rank approximation of S_{k}(f,g)
        [fx_lr, gx_lr, hx_lr, ux_lr, vx_lr, wx_lr, lambda_lr, theta_lr] =...
            SNTLN_3Polys(fx, gx, hx, lambda, theta, k, idx_col);
        
        % Get f(\omega)
        lambda_fw   = lambda_lr .* GetWithThetas(fx, theta_lr);
        mu_gw       = mu_lr     .* GetWithThetas(gx, theta_lr);
        rho_hw      = rho_lr    .* GetWithThetas(hx, theta_lr)
        
        % Plot the Singular values of the Sylvester subresultants S_{k}
        Plot_LowRank_SingularValues(fx, gx, fx_lr, gx_lr,...
            lambda_fw, mu_gw, lambda_fw, mu_gw,k)
        
        
        
        
        
        
    case 'None'
        
        
        % Get f(\omega) and g(\omega)
        lambda_fw   = lambda    .* GetWithThetas(fx, theta);
        mu_gw       = mu        .* GetWithThetas(gx,  theta);
        rho_hw      = rho       .* GetWithThetas(hx,  theta);
        
        
        
        % Get quotient polynomials 
        %   \tilde{u}_{k}(\omega), 
        %   \tilde{v}_{k}(\omega) and 
        %   \tilde{w}_{k}(\omega)
        [uw, vw, ww] = GetQuotients_3Polys(lambda_fw, mu_gw, rho_hw, k);
        
        
        % f(x), g(x) and h(x) are unchanged by this function
        fx_lr = fx;       
        gx_lr = gx;
        hx_lr = hx;
        
        % \lambda, \mu, \rho and \theta are unchanged by this function
        lambda_lr   = lambda;
        mu_lr       = mu;
        rho_lr      = rho;
        theta_lr    = theta;
        
        % Get polynomials u(x), v(x) and w(x) from u(\omega), v(\omega)
        % and w(\omega)
        vx_lr = GetWithoutThetas(vw, theta);
        ux_lr = GetWithoutThetas(uw, theta);
        wx_lr = GetWithoutThetas(ww, theta);
        
        SETTINGS.LOW_RANK_APPROX_REQ_ITE = 0;
        
    otherwise
        
        error('SETTINGS.LOW_RANK_APPROXIMATION_METHOD must be valid')
        
end

end


function Plot_LowRank_SingularValues(fx, gx, hx,...
    fx_lr, gx_lr, hx_lr, ...
    fw, gw, hw, ...
    fw_lr, gw_lr, hw_lr, ...
    k)
%
% % Inputs
%
% fx : (Vector) The coefficients of the input polynomial f(x)
%
% gx : (Vector) The coefficients of the input polynomial g(x)
%
% hx : (Vector) The coefficients of the input polynomial h(x)
%
% fx_lr : (Vector) The Coefficients of polynomial f(x) + \delta f(x) from low rank
% approximation method.
%
% gx_lr : (Vector) The coefficients of polynomial g(x) + \delta g(x) from 
% a low rank approximation method
%
% hx_lr : (Vector) The coefficients of polynomial h(x) + \delta h(x) from 
% a low rank approximation method
%
% fw : (Vector) The coefficients of the polynomial f(\omega)
%
% gw : (Vector) The coefficients of the polynomial g(\omega)
%
% hw : (Vector) The coefficients of the polynomial h(\omega)
%
% fw_lr : (Vector) The coefficients of the polynomial f(\omega) LR
%
% gw_lr : (Vector) The coefficients of the polynomial g(\omega) LR
%
% hw_lr : (Vector) The coefficients of the polynomial h(\omega) LR


global SETTINGS

if(SETTINGS.PLOT_GRAPHS)
    
    
    % Build Sylvester subresultant matrices S_{k} for each of the four
    % pairs of polynomials.
    
    S1 = BuildDTQ_3Polys(fx, gx, hx, k);
    S2 = BuildDTQ_3Polys(fx_lr, gx_lr, hx_lr, k);
    S3 = BuildDTQ_3Polys(fw, gw, hw, k);
    S4 = BuildDTQ_3Polys(fw_lr, gw_lr, hw_lr, k);
    
    % Get singular values for each of the 4 Sylvester subresultant
    % matrices S_{k}
    vSingularValues1 = svd(S1);
    vSingularValues2 = svd(S2);
    vSingularValues3 = svd(S3);
    vSingularValues4 = svd(S4);
    
    % Plot Singular values
    figure_name = sprintf([mfilename ' : Singular Values']);
    figure('name',figure_name)
    
    plot(log10(vSingularValues1),'-s','DisplayName', ...
        '$f(x),g(x),h(x)$')
    hold on
    plot(log10(vSingularValues2),'-s','DisplayName', ...
        '$f(x),g(x),h(x) LR$')
    plot(log10(vSingularValues3),'-s','DisplayName', ...
        '$\tilde{f}(\omega),\tilde{g}(\omega), \tilde{h}(\omega)$')
    plot(log10(vSingularValues4),'-s','DisplayName', ...
        '$f(\omega), \tilde{g}(\omega), \tilde{h}(\omega)$')
    
    
    l = legend(gca,'show');
    set(l,'Interpreter','latex')
    
    hold off
    
end

end
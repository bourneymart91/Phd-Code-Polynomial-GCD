function [fx_lr, gx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = ...
    LowRankApproximation_2Polys(fx, gx, GM_fx, GM_gx, alpha, theta, k, idx_col)
% Get the low rank approximation of the Sylvester subresultant matrix
% S_{k}(f,g)
%
% Inputs.
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x)
%
% alpha :  (Float) Optimal value of \alpha
%
% theta : (Float) Optimal value of \theta
%
% k : (Int) Degree of GCD of f(x) and g(x)
%
% idx_col : (Int) Index of optimal column for removal from S_{k}(f,g)
%
% % Outputs
%
% fx_lr : (Vector) Output coefficients of the polynomial f(x)
%
% gx_lr : (Vector) Output coefficients of the polynomial g(x)
%
% ux_lr : (Vector) Output coefficients of the polynomial u(x)
%
% vx_lr : (Vector) Output coefficients of the polynomial v(x)
%
% alpha_lr : (Float) Optimal value of \alpha
%
% theta_lr : (Float) Optimal value of \theta

global SETTINGS

switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    
  
        
        
        
    
    case 'Standard STLN'
        
        % Get preprocessed polynomials f(\omega) and g(\omega)
        
        fx_n = fx ./ GM_fx;
        gx_n = gx ./ GM_gx;
        
        fw = GetWithThetas(fx_n, theta);
        gw = GetWithThetas(gx_n, theta);
        a_gw = alpha .* gw;
        
        % Performe STLN to get low rank approximation of S(f,g)
        [fw_lr, a_gw_lr, uw_lr, vw_lr] = STLN(fw, a_gw, k, idx_col);
        
        % Get f(x) and g(x) from low rank approximation.
        fx_lr = GetWithoutThetas(fw_lr, theta);
        gx_lr = GetWithoutThetas(a_gw_lr, theta) ./ alpha;
        
        % Get u(x) and v(x) from low rank approximation
        ux_lr = GetWithoutThetas(uw_lr, theta);
        vx_lr = GetWithoutThetas(vw_lr, theta);
        
        % \alpha and \theta are same as input, as they are unchanged by
        % STLN
        alpha_lr = alpha;
        theta_lr = theta;
        
        %
        Plot_LowRank_SingularValues(fx, gx, fx_lr, gx_lr,...
            fw, alpha.*gw, fw_lr, a_gw_lr,k);
        
        
        %PlotCoefficients({fw, fw_lr},{'fw','fw LR'})
        %PlotCoefficients({a_gw, a_gw_lr},{'gw','gw LR'})
        
        %PlotCoefficients({fx, fx_lr},{'fx','fx LR'})
        %PlotCoefficients({gx, gx_lr},{'gx','gx LR'})
        
        %PlotCoefficients({fx, fx_lr, fw, fw_lr},{'fx','fx LR','fw','fw LR'})
        %PlotCoefficients({gx, gx_lr, alpha.*gw, a_gw_lr},{'gx','gx LR','a_gw','a_gw LR'})
        
        
        % Residual for fx gx ux vx could be large when alpha is not
        % included. Alpha must be included in computing the residual. An
        % alternative method would be to consider the 2norm error between
        % fv and gu
        
        res1 = GetResidual(fx, alpha.* gx, ux_lr, vx_lr, k);
        res2 = GetResidual(fx_lr, alpha .*gx_lr, ux_lr, vx_lr, k);
        res3 = GetResidual(fw, a_gw, uw_lr, vw_lr, k);
        res4 = GetResidual(fw_lr, a_gw_lr, uw_lr, vw_lr, k);
        
        err1 = GetError(fx, alpha.*gx, ux_lr, vx_lr, k);
        err2 = GetError(fx_lr, alpha.*gx_lr, ux_lr, vx_lr, k);
        err3 = GetError(fw, a_gw, uw_lr, vw_lr, k);
        err4 = GetError(fw_lr, a_gw_lr, uw_lr, vw_lr, k);
        
        fprintf('Residual f(x), g(x), u(x), v(x)    : %g \n', res1)
        fprintf('Residual f(x), g(x), u(x), v(x) LR : %g \n', res2)
        fprintf('Residual f(w), g(w), u(w), v(w)    : %g \n', res3)
        fprintf('Residual f(w), g(w), u(w), v(w) LR : %g \n', res4)
        
        fprintf('Error 1 : %g \n', err1)
        fprintf('Error 2 : %g \n', err2)
        fprintf('Error 3 : %g \n', err3)
        fprintf('Error 4 : %g \n', err4)
        
    case 'Standard SNTLN' % Structured Non-Linear Total Least Norm
        
        % Perform SNTLN to get low rank approximation of S_{k}(f,g)
        fx_n = fx ./ GM_fx;
        gx_n = gx ./ GM_gx;
        
        [fx_lr, gx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = SNTLN(fx_n, gx_n, alpha, theta, k, idx_col);
        
        % Get f(\omega)
        fw = GetWithThetas(fx, theta_lr);
        
        % Get g(\omega)
        a_gw = alpha_lr.* GetWithThetas(gx, theta_lr);
        
        % Get f(\omega) from f(x) from low rank approximation
        fw_lr = GetWithThetas(fx_lr, theta_lr);
        
        % Get g(\omega) from g(x) from low rank approximation
        a_gw_lr = alpha_lr .* GetWithThetas(gx_lr, theta_lr);
        
        
        uw_lr = GetWithThetas(ux_lr, theta_lr);
        vw_lr = GetWithThetas(vx_lr, theta_lr);
        
        
        % Plot the Singular values of the Sylvester subresultants S_{k}
        Plot_LowRank_SingularValues(fx, gx, fx_lr, gx_lr,...
            fw, a_gw, fw_lr, a_gw_lr,k)
        
        
        %PlotCoefficients({fw, fw_lr},{'$f(\omega)$','$f(\omega) LR$'},{'-r','-b'})
        %PlotCoefficients({a_gw, a_gw_lr},{'$g(\omega)$','$g(\omega) LR$'},{'-r','-b'})
        
        %PlotCoefficients({fx, fx_lr},{'fx','fx LR'})
        %PlotCoefficients({gx, gx_lr},{'gx','gx LR'})
        
        %PlotCoefficients({fx, fx_lr, fw, fw_lr},{'fx','fx LR','fw','fw LR'})
        %PlotCoefficients({gx, gx_lr, a_gw, a_gw_lr},{'gx','gx LR','a_gw','a_gw LR'})
        
        
        % Residual for fx gx ux vx could be large when alpha is not
        % included. Alpha must be included in computing the residual. An
        % alternative method would be to consider the 2norm error between
        % fv and gu
        
        res1 = GetResidual(fx, alpha.* gx, ux_lr, vx_lr, k);
        res2 = GetResidual(fx_lr, alpha_lr .*gx_lr, ux_lr, vx_lr, k);
        res3 = GetResidual(fw, a_gw, uw_lr, vw_lr, k);
        res4 = GetResidual(fw_lr, a_gw_lr, uw_lr, vw_lr, k);
        
        err1 = GetError(fx, alpha.*gx, ux_lr, vx_lr, k);
        err2 = GetError(fx_lr, alpha_lr.*gx_lr, ux_lr, vx_lr, k);
        err3 = GetError(fw, a_gw, uw_lr, vw_lr, k);
        err4 = GetError(fw_lr, a_gw_lr, uw_lr, vw_lr, k);
        
        fprintf('Residual f(x), g(x), u(x), v(x)    : %g \n', res1)
        fprintf('Residual f(x), g(x), u(x), v(x) LR : %g \n', res2)
        fprintf('Residual f(w), g(w), u(w), v(w)    : %g \n', res3)
        fprintf('Residual f(w), g(w), u(w), v(w) LR : %g \n', res4)
        
        fprintf('Error 1 : %g \n', err1)
        fprintf('Error 2 : %g \n', err2)
        fprintf('Error 3 : %g \n', err3)
        fprintf('Error 4 : %g \n', err4)
        
        
    case 'STLN Root Specific'
        
        % Exclude the geometric mean until this has been thought through in
        % STLNRootSpecific
        
        fx_n = fx ./ GM_fx;
        gx_n = gx ./ GM_gx;
        
        % Get preprocessed polynomials f(\omega) and g(\omega)
        fw = GetWithThetas(fx_n, theta);
        gw = GetWithThetas(gx_n, theta);
        a_gw = alpha.* gw;
        
        % Performe STLN to get low rank approximation of S(f,g)
        [fw_lr, a_gw_lr, uw_lr, vw_lr] = STLNRootSpecific(fx, gx, GM_fx, GM_gx, alpha, theta, k, idx_col);
        
        % Get f(x) and g(x) from low rank approximation.
        fx_lr = GetWithoutThetas(fw_lr, theta);
        gx_lr = GetWithoutThetas(a_gw_lr, theta) ./ alpha;
        
        % Get u(x) and v(x) from low rank approximation
        ux_lr = GetWithoutThetas(uw_lr, theta);
        vx_lr = GetWithoutThetas(vw_lr, theta);
        
        % \alpha and \theta are same as input, as they are unchanged by
        % STLN
        alpha_lr = alpha;
        theta_lr = theta;
        
        %
        Plot_LowRank_SingularValues(fx, gx, fx_lr, gx_lr,...
            fw, alpha.*gw, fw_lr, a_gw_lr,k);
        
        
        %PlotCoefficients({fw, fw_lr},{'fw','fw LR'})
        %PlotCoefficients({a_gw, a_gw_lr},{'gw','gw LR'})
        
        %PlotCoefficients({fx, fx_lr},{'fx','fx LR'})
        %PlotCoefficients({gx, gx_lr},{'gx','gx LR'})
        
        %PlotCoefficients({fx, fx_lr, fw, fw_lr},{'fx','fx LR','fw','fw LR'})
        %PlotCoefficients({gx, gx_lr, alpha.*gw, a_gw_lr},{'gx','gx LR','a_gw','a_gw LR'})
        
        
        % Residual for fx gx ux vx could be large when alpha is not
        % included. Alpha must be included in computing the residual. An
        % alternative method would be to consider the 2norm error between
        % fv and gu
        
        res1 = GetResidual(fx, alpha.* gx, ux_lr, vx_lr, k);
        res2 = GetResidual(fx_lr, alpha .*gx_lr, ux_lr, vx_lr, k);
        res3 = GetResidual(fw, a_gw, uw_lr, vw_lr, k);
        res4 = GetResidual(fw_lr, a_gw_lr, uw_lr, vw_lr, k);
        
        err1 = GetError(fx, alpha.*gx, ux_lr, vx_lr, k);
        err2 = GetError(fx_lr, alpha.*gx_lr, ux_lr, vx_lr, k);
        err3 = GetError(fw, a_gw, uw_lr, vw_lr, k);
        err4 = GetError(fw_lr, a_gw_lr, uw_lr, vw_lr, k);
        
        fprintf('Residual f(x), g(x), u(x), v(x)    : %g \n', res1)
        fprintf('Residual f(x), g(x), u(x), v(x) LR : %g \n', res2)
        fprintf('Residual f(w), g(w), u(w), v(w)    : %g \n', res3)
        fprintf('Residual f(w), g(w), u(w), v(w) LR : %g \n', res4)
        
        fprintf('Error 1 : %g \n', err1)
        fprintf('Error 2 : %g \n', err2)
        fprintf('Error 3 : %g \n', err3)
        fprintf('Error 4 : %g \n', err4)
        
    case 'SNTLN Root Specific'
        
        error('Code Not Complete')
        
    case 'None'
        
        
        % Get f(\omega) and g(\omega)
        fw = GetWithThetas(fx, theta);
        gw = GetWithThetas(gx, theta);
        
        % Get quotient polynomials u(\omega) and v(\omega)
        [uw_lr, vw_lr] = GetQuotients_2Polys(fw, alpha.*gw, k);
        
        % f(x) and g(x) are unchanged from input
        fx_lr = fx;
        gx_lr = gx;
        
        % \alpha and \theta are unchanged from input
        alpha_lr = alpha;
        theta_lr = theta;
        
        % Get polynomials u(x) and v(x) from u(\omega) and v(\omega)
        vx_lr = GetWithoutThetas(vw_lr, theta);
        ux_lr = GetWithoutThetas(uw_lr, theta);
        SETTINGS.LOW_RANK_APPROX_REQ_ITE = 0;
    otherwise
        
        error('SETTINGS.LOW_RANK_APPROXIMATION_METHOD must be valid')
        
end

end


function Plot_LowRank_SingularValues(fx, gx, fx_lr, gx_lr, fw, a_gw, fw_lr, a_gw_lr,k )
%
% % Inputs
%
% fx : (Vector) Coefficients of input polynomial f(x)
%
% gx : (Vector) Coefficients of input polynomial g(x)
%
% fx_lr : (Vector) Coefficients of polynomial f(x) + \delta f(x) from low rank
% approximation method.
%
% gx_lr : (Vector) Coefficients of polynomial g(x) + \delta g(x) from low rank
% approximation method
%
% fw : (Vector) Coefficients of polynomial f(\omega)
%
% gw : (Vector) Coefficients of polynomial g(\omega)
%
% fw_lr : (Vector)
%
% gw_lr : (Vector)

global SETTINGS

if(SETTINGS.PLOT_GRAPHS_LOW_RANK_APPROXIMATION)
    
    
    bool_normalise = true;
    bool_log = true;
    
    % Build Sylvester subresultant matrices S_{k} for each of the four
    % pairs of polynomials.
    
    % Unprocessed
    arrSk{1} = BuildDTQ(fx, gx, 1);
    
    % Unprocessed after low rank approximation
    arrSk{2} = BuildDTQ(fx_lr, gx_lr, 1);
    
    % Preprocessed
    arrSk{3} = BuildDTQ(fw, a_gw, 1);
    
    % Preprocessed after low rank approximation
    arrSk{4} = BuildDTQ(fw_lr, a_gw_lr, 1);
    
    % Get singular values for each of the 4 Sylvester subresultant
    % matrices S_{k}
    
    nSylvesterMatrices = 4;
    
    arrLabel = {...
        'S(f(x), g(x))',...
        'S(f(x),g(x)) LR',...
        'S(f(\omega), g(\omega))',...
        'S(f(\omega), g(\omega)) LR'};
    
    arrVSingularValues = cell(nSylvesterMatrices,1);
    
    
    for i = 1:1:nSylvesterMatrices
        
        temp_vec = svd(arrSk{i});
        
        
        
        if bool_normalise == true
            temp_vec = temp_vec./ temp_vec(1);
        end
        
        arrVSingularValues{i} = temp_vec;
        
    end
    
    
    
    % Plot Singular values
    figure_name = sprintf([mfilename ' : Singular Values']);
    figure('name',figure_name)
    hold on

    for i = 1 : 1 : nSylvesterMatrices
    
        vSingularValues = arrVSingularValues{i};
        strLabel = arrLabel{i};
 
        
        if bool_log == true
            vSingularValues = log10(vSingularValues);
        end
            
        plot(vSingularValues, '-s', 'DisplayName', strLabel,'LineWidth',2)
        
    end

    legend(gca,'show');
    hold off
    
    
end
end


function residual = GetResidual(fx, gx, ux, vx, k)


Sk = BuildSubresultant_2Polys(fx,gx,k);
vec_x = [vx;-1.* ux];

residual = Sk*vec_x;

residual = norm(residual);

end

function err = GetError(fx,gx,ux,vx,k)

m = GetDegree(fx);
n = GetDegree(gx);

Cf = BuildSubresultant_Partition_2Polys(fx, n - k);
Cg = BuildSubresultant_Partition_2Polys(gx, m - k);

fv = Cf * vx;
gu = Cg * ux;


err = norm(fv - gu) ./ norm(fv);

end



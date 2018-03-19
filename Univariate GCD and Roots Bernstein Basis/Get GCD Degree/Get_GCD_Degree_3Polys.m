function [t, lambda, mu, rho, theta, GM_fx, GM_gx, GM_hx] = ...
    Get_GCD_Degree_3Polys(fx, gx, hx, limits_t, rank_range)
% GetGCD_Degree_2Polys(fx,gx)
%
% Get degree t of the AGCD d(x) of input polynomials f(x) and g(x)
%
%
% % Inputs.
%
% fx : (Vector) Vector of the coefficients of the polynomial f(x)
%
% gx : (Vector) Vector of the coefficients of the polynomail g(x)
%
% hx : (Vector) Vector of the coefficients of the polynomial h(x)
%
% limits_t : [Int Int] Set the upper and lower bound of the degree of the
% GCD of polynomials f(x) and g(x). Usually used when using o_roots() where
% prior information is known about the upper and lower bound due to number
% of distinct roots.
%
% rank_range : [Float Float] : Used as upper and lower threshold for
% determining whether a matrix is singular or non-singular.
%
%
% % Outputs.
%
% t : (Int) Degree of GCD of f(x) and g(x)
%
% lambda : (Float) Optimal value of \lambda_{t} in preprocessing the
% coefficients of f(x) in the t-th subresultant matrix S_{t}(f,g,h)
%
% mu : (Float) Optimal value of \mu_{t} in preprocessing the coefficients
% of the polynomial g(x) in the t-th subresultant matrix S_{t}(f,g,h)
%
% rho : (Float) Optimal value of \rho_{t} in preprocessing the coefficients
% of the polynomial h(x) in the t-th subresultant matrix S_{t}(f,g,h)
%
% theta : (Float) optimal value of \theta in the t-th subresultant matrix
%
% GM_fx : (Float) Geometric mean of the entries of f(x) in the t-th
% subresultant matrix
%
% GM_gx : (Float) Geometric mean of the entries of g(x) in the t-th
% subresultant matrix
%
% GM_hx : (Float) Geometric mean of the entries of h(x) in the t-th 
% subresultant matrix 


global SETTINGS




% if not 5 input arguments, then error
if (nargin ~= 5)
    error('Not enough input arguments');
end




% Get the degree of the polynomails f(x), g(x) and h(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

% Set upper and lower limit for computing the degree of the GCD. Note may
% be best to set to degree limits or may be best to set to 1 & min(m,n)
lowerLimit_k = 1;
upperLimit_k = min([m, n, o]);
limits_k = [lowerLimit_k upperLimit_k];



% Get the number of subresultants which must be constructed - Note that we
% always construct all subresultant matrices in the general GCD finding 
% problem, however, limits can be defined in the GCD problem in the
% polynomial square free factorisation algorithm.
nSubresultants = upperLimit_k - lowerLimit_k + 1 ;




% %
% Initialisation stage

% Initialise vectors to store all optimal alphas and theta, and each
% geometric mean for f(x), g(x) and h(x) in each S_{k} for k = 1,...,min(m,n)
vLambda = zeros(nSubresultants, 1);
vMu = zeros(nSubresultants, 1);
vRho = zeros(nSubresultants, 1);
vTheta = zeros(nSubresultants, 1);
vGM_fx = zeros(nSubresultants, 1);
vGM_gx = zeros(nSubresultants, 1);
vGM_hx = zeros(nSubresultants, 1);

% Initialise arrays to store Sylvester matrices
arr_Sk = cell(nSubresultants, 1);
arr_R1 = cell(nSubresultants, 1);

% Initialise vectors
vMaxDiagR1 = zeros(nSubresultants,1);
vMinDiagR1 = zeros(nSubresultants,1);
vMaxRowNormR1 = zeros(nSubresultants,1);
vMinRowNormR1 = zeros(nSubresultants,1);

% Initialise a vector to store minimal residuals obtained by QR
% decomposition of each subresultant S_{k} for k=1,...,min(m,n)
vMinimumResidual_QR  = zeros(nSubresultants,1);
vMinimumResidual_SVD = zeros(nSubresultants,1);

% Initialise a vector to store the minimum singular values \sigma_{i} of each
% subresultant S_{i}(f,g), where i is between the upper and lower bound.
vMinimumSingularValues = zeros(nSubresultants,1);

vCondition = zeros(nSubresultants,1);


% For each subresultant $S_{k}$
for k = lowerLimit_k : 1 : upperLimit_k
    
    i = k - lowerLimit_k + 1;
    
    
    
    switch SETTINGS.SYLVESTER_EQUATIONS
        
        case '2'
            [GM_fx, GM_gx, GM_hx, lambda, mu, rho, theta ] = ...
                Preprocess_3Polys_2Eqns(fx, gx, hx, k);
            
        case '3'
            [GM_fx, GM_gx, GM_hx, lambda, mu, rho, theta ] = ...
                Preprocess_3Polys_3Eqns(fx, gx, hx, k);
        otherwise
            error('err')
            
    end
    
    % Store geometric means of entries of f(x), g(x) and h(x) found in the
    % i-th subresultant matrix
    vGM_fx(i) = GM_fx;
    vGM_gx(i) = GM_gx;
    vGM_hx(i) = GM_hx;
    
    % Store optimal values of \lambda, \mu, \rho and \theta in the i-th
    % subresultant matrix
    vLambda(i) = lambda;
    vMu(i) = mu;
    vRho(i) = rho;
    vTheta(i) = theta;
    
    % Divide f(x) and g(x) by geometric means
    fx_n = fx ./ vGM_fx(i);
    gx_n = gx ./ vGM_gx(i);
    hx_n = hx ./ vGM_hx(i);
    
    % Construct the i-th preprocessed subresultant matrix
    % S_{k}(f(\theta),g(\theta))
    alpha_fw = lambda .* GetWithThetas(fx_n, vTheta(i));
    beta_gw = mu .* GetWithThetas(gx_n, vTheta(i));
    gamma_hw = rho .* GetWithThetas(hx_n, vTheta(i));
    
    % Build the i-th subresutlant matrix
    arr_Sk{i} = BuildSubresultant_3Polys(alpha_fw, beta_gw, gamma_hw, k);
    
    
    
    
    if k == lowerLimit_k
        
        if (SETTINGS.PLOT_GRAPHS_PREPROCESSING == true)
            
            figure_name = sprintf('Heat Map : %s : %s', ...
                SETTINGS.SYLVESTER_MATRIX_VARIANT, ...
                SETTINGS.EX_NUM);
            
            figure('Name',figure_name)
            hold on
            colormap('hot');
            data = log10(abs(arr_Sk{i}));
            imagesc(flipud(data));
            colorbar;
            hold off
            
            
            PlotCoefficients(...
                {...
                fx, alpha_fw,...
                gx, beta_gw,...
                hx, gamma_hw...
                },...
                {...
                '$f(x)$', '$\lambda \tilde{f}(\omega)$',...
                '$g(x)$', '$\mu \tilde{g}(\omega)$',...
                '$h(x)$', '$\rho \tilde{h}(\omega)$'...
                },...
                {'--','-s','--','-s','--','-s'}...
                );
            
        end
    end
    
    
    
    vCondition(i) = cond(arr_Sk{i});
    
    % Get the matrix R1 from the QR Decomposition of S
    arr_R1{i} = GetR1(arr_Sk{i});
    
    
    
    
end % End of for







if (SETTINGS.PLOT_GRAPHS_PREPROCESSING == true)
    
    % % Plot Lambda, Mu, rho and theta
    figure_name = strcat('Scaling : ',  SETTINGS.SCALING_METHOD);
    figure('Name',figure_name);
    hold on
    
    plot(log10(vLambda), '-s', 'LineWidth', 2, 'DisplayName', '\lambda')
    plot(log10(vMu), '-o', 'LineWidth', 2, 'DisplayName', '\mu')
    plot(log10(vRho), '-*', 'LineWidth', 2, 'DisplayName', '\rho')
    plot(log10(vTheta), 'LineWidth', 2, 'DisplayName', '\theta')
    
    ylabel('$\log_{10} \left( \Re \right)$', 'Interpreter','latex')
    xlabel('$k$', 'Interpreter','latex')
    l = legend(gca,'show');
    set(l,'location','southwest')
    set(l,'FontSize',20)
    hold off
    
    
    % Plot Geometric Means
    figure('Name','Geometric Means')
    hold on
    
    plot(log10(vGM_fx), '-s', 'DisplayName','$GM f(x)$')
    plot(log10(vGM_gx), '-o', 'DisplayName','$GM g(x)$')
    plot(log10(vGM_hx), '-', 'DisplayName','$GM h(x)$')
    l = legend(gca,'show');
    set(l,'Interpreter', 'latex')
    hold off
    
    % Plot Condition number
    figure('Name','Condition Numbers')
    hold on
    plot(log10(vCondition));
    hold off
    
end



% %
% %
% %
% Choose a metric to determine the degree of the GCD.
% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals


fprintf(sprintf('Metric used to compute degree of GCD : %s \n', SETTINGS.RANK_REVEALING_METRIC));

switch SETTINGS.RANK_REVEALING_METRIC
    case 'Row Norms'
        
        for i = 1:1:nSubresultants
            
            % Get Norms of each row in the matrix R1
            vR1_RowNorms = sqrt(sum(arr_R1{i}.^2,2))./norm(R1);
            
            % Get maximum and minimum row norms of rows of R1.
            vMaxRowNormR1(i) = max(vR1_RowNorms);
            vMinRowNormR1(i) = min(vR1_RowNorms);
            
        end
        
        vMetric = log10(vMaxRowNormR1./vMinRowNormR1);
        
        plotMaxMinRowSum(vMaxRowNormR1, vMinRowNormR1, limits_k, limits_t);
        
    case 'Row Diagonals'
        
        for i = 1:1:nSubresultants
            
            % Get Norms of diagonal eleements of R1
            vDiagsR1 = diag(arr_R1{i});
            
            % Get maximum and minimum row diagonals of R1
            vMaxDiagR1(i) = max(vDiagsR1);
            vMinDiagR1(i) = min(vDiags1);
            
        end
        
        plotMaxMinRowDiagonals(vMaxDiagR1, vMinDiagR1, limits_k, limits_t);
        
        vMetric = log10(vMaxDiagR1./vMinDiagR1);
        
    case 'Minimum Singular Values'
        
        arr_SingularValues = cell(nSubresultants,1);
        
        for i = 1 : 1 : nSubresultants
            
            % Get singular values of S_{k}
            arr_SingularValues{i} = svd(arr_Sk{i});
            
            % Get the minimal singular value from S_{k}
            vMinimumSingularValues(i) = min(arr_SingularValues{i});
            
        end
        
        vMetric = log10(vMinimumSingularValues);
        
        plotSingularValues(arr_SingularValues, limits_k, limits_t);
        %plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t, rank_range);
        
    case 'Normalised Minimum Singular Values'
        
        arr_NormalisedSingularValues = cell(nSubresultants,1);
        
        for i = 1 : 1 : nSubresultants
            
            % Get singular values of S_{k}
            
            vSingularValues = svd(arr_Sk{i});
            vNormalisedSingularValues = vSingularValues./vSingularValues(1);
            
            arr_NormalisedSingularValues{i} = vNormalisedSingularValues;
            
            % Get the minimal singular value from S_{k}
            vNormalisedMinimumSingularValues(i) = min(arr_NormalisedSingularValues{i});
            
        end
        
        vMetric = log10(vNormalisedMinimumSingularValues);
        
        plotSingularValues(arr_NormalisedSingularValues, limits_k, limits_t);
        plotMinimumSingularValues(vNormalisedMinimumSingularValues, limits_k, limits_t, rank_range);
        
    case 'Residuals'
        
        % Get the minimal residuals for each subresultant S_{k}.
        vMinimumResidual_QR(i) = GetMinimalDistance(arr_Sk{i},'QR');
        vMinimumResidual_SVD(i) = GetMinimalDistance(arr_Sk{i},'SVD');
        
        error('Code not developed')
        
    otherwise
        error('err')
end








% % Analysis of Minimum Singular values

if (upperLimit_k == lowerLimit_k)
    
    % Set \alpha and \theta
    lambda = vLambda(1);
    mu = vMu(1);
    theta = vTheta(1);
    
    % Set Geometric means
    GM_fx = vGM_fx(1);
    GM_gx = vGM_gx(1);
    GM_hx = vGM_hx(1);
    
    % Set degree of GCD
    t = upperLimit_k;
    
    return;
end










% If only one subresultant exists, use an alternative method.
if (upperLimit_k - lowerLimit_k == 0 )
    
    % Set the degree of the GCD
    t = Get_GCD_Degree_OneSubresultant(vSingularValues);
    
    % Set \alpha and \theta
    lambda = vLambda(1);
    mu = vMu(1);
    theta = vTheta(1);
    
    % Set geometric means
    GM_fx = vGM_fx(1);
    GM_gx = vGM_gx(1);
    GM_hx = vGM_hx(1);
    
    return;
    
    
    
else
    
    [t] = Get_GCD_Degree_MultipleSubresultants_2Polys(vMetric, limits_k, limits_t, rank_range);
    
    
    % % Graph Plotting
    %PlotGraphs_GCDDegree()
    
    % %
    % Outputs
    
    % Output all subresultants, all optimal alphas, all optimal thetas and all
    % geometric means for each subresultant S_{k} where k = 1,...,min(m,n)
    
    lambda = vLambda(t);
    theta = vTheta(t);
    
    GM_fx = vGM_fx(t);
    GM_gx = vGM_gx(t);
    GM_hx = vGM_hx(t);
    
end




end


function data = AddToResults(data,vector,k)
% Given the vector of data, it is necessary to include in a 3 columned
% matrix, along with the corresponding value k to indicate that the value
% came from the k-th subresultant, and an index 1:1:nEntries.


% Get Number of entries in vector
nEntries = size(vector,1);

% Get a vector of k
v_k = k.* ones(nEntries,1);

% Get a vector of 1:1:n
vec_n = 1:1:nEntries;

% Form a triple of [ks, the value of QR_RowNorm, and the index of the value of
% the row of R1 corresponding to QR_RowNorm].

temp_data = [v_k vector vec_n'];
data = [data; temp_data];

end





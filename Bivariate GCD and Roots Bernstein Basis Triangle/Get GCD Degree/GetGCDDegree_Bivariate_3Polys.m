function [t, GM_fx, GM_gx, GM_hx, lambda, mu, rho, th1, th2, rank_range] = ...
    GetGCDDegree_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, limits_t, rank_range)
% GetGCDDegree_Bivariate_3Polys(fxy, gxy, hxy, m, n, o)
%
% % Inputs.
%
% fxy : (Vector) Coefficients of polynomial f(x,y)
%
% gxy : (Vector) Coefficients of polynomial g(x,y)
%
% hxy : (Vector) Coefficients of polynomial h(x,y)
%
% m : (Int) Degree of polynomial f(x,y)
%
% n : (Int) Degree of polynomial g(x,y)
%
% o : (Int) Degree of polynomial h(x,y)
%
% limits_t : [(Int) (Int)]
%
% rank_range : [(Float) (Float) ]
%
% % Outputs
%
% t : (Int) Degree of GCD d(x,y)
%
% GM_fx1 : (Float) Geometric mean of coefficients of f(x) in the first
% partition of the subresultant matrix.
%
% GM_fx2 : (Float) Geometric mean of coefficients of f(x) in the third
% partition of the subresultant matrix.
%
% GM_gx : (Float) Geometric mean of coefficients of g(x)
%
% GM_hx : (Float) Geometric mean of coefficients of h(x)
%
% alpha : (Float) Optimal value of \alpha
%
% beta : (Float) Optimal value of \beta
%
% gamma : (Float) Optimal value of \gamma
%
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}

global SETTINGS

% Get range of Sylvester subresultant matrices to be constructed
lowerLimit_k = 1;
upperLimit_k = min([m,n,o]);
limits_k = [lowerLimit_k upperLimit_k];

rank_range_low = rank_range(1);
rank_range_high = rank_range(2);

% Get number of Sylvester Subresultants
nSubresultants = upperLimit_k - lowerLimit_k + 1;

arr_SingularValues = cell(nSubresultants, 1);
arr_NormalisedSingularValues = cell(nSubresultants, 1);
arr_R1 = cell(nSubresultants, 1);

% Initialise vectors to store geometric means
vGM_fx = zeros(nSubresultants, 1);
vGM_gx = zeros(nSubresultants, 1);
vGM_hx = zeros(nSubresultants, 1);

% Initialise vectors to store alpha, beta, theta_{1} and theta_{2}
vAlpha = zeros(nSubresultants, 1);
vBeta = zeros(nSubresultants, 1);
vGamma = zeros(nSubresultants, 1);

vTh1 = zeros(nSubresultants, 1);
vTh2 = zeros(nSubresultants, 1);

for i = 1 : 1 : nSubresultants
    
    k = lowerLimit_k + (i-1);
    
    % %
    % Geometric Mean
    
    
    switch SETTINGS.SYLVESTER_MATRIX_3POLY_N_EQUATIONS
        
        case '2'
            % Get Geometric mean of f(x,y) in C_{n-k}(f), g(x,y) and h(x,y)
            vGM_fx(i) = GetMean_2Partitions(fxy, m, n, o, k );
            vGM_gx(i) = GetMean(gxy, n, m - k);
            vGM_hx(i) = GetMean(hxy, o, n - k);
            
        case '3'
            vGM_fx(i) = GetMean_2Partitions(fxy, m, n, o, k );
            vGM_gx(i) = GetMean_2Partitions(gxy, n, m, o, k);
            vGM_hx(i) = GetMean_2Partitions(hxy, o, n, m, k);
            
        otherwise
            error('err')
    end
    
    % Divide entries of f(x,y) and g(x,y) by geometric mean
    fxy_n = fxy ./ vGM_fx(i);
    gxy_n = gxy ./ vGM_gx(i);
    hxy_n = hxy ./ vGM_hx(i);
    
    
    % %
    % Preprocess
    [lambda, mu, rho, th1, th2] = Preprocess_3Polys(fxy_n, gxy_n, hxy_n, m, n, o, k);
    
    
    
    
    
    vAlpha(i) = lambda;
    vBeta(i) = mu;
    vGamma(i) = rho;
    vTh1(i) = th1;
    vTh2(i) = th2;
    
    % Get f(x,y) and g(x,y) with thetas to get f(w_{1},w_{2}) and g(w_{1},w_{2})
    fww = GetWithThetas(fxy_n, m, th1, th2);
    gww = GetWithThetas(gxy_n, n, th1, th2);
    hww = GetWithThetas(hxy_n, o, th1, th2);
    
    if i == 1
        PlotCoefficients(...
            {fxy, lambda .* fww,...
            gxy, mu.*gww, ...
            hxy, rho.*hww},...
            {...
            '$f(x,y)$', '$\lambda \tilde{f}(\omega_{1},\omega_{2})$', ...
            '$g(x,y)$', '$\mu     \tilde{g}(\omega_{1}, \omega_{2})$', ...
            '$h(x,y)$', '$\rho    \tilde{h}(\omega_{1}, \omega_{2})$'...
            });
        %PlotCoefficients({gxy, mu .* gww}, {'g(x,y)', '\beta g(\omega, \omega)'});
        %PlotCoefficients({hxy, rho .* hww}, {'h(x,y)', '\gamma h(\omega,\omega)'});
    end
    
    
    % Build the Sylvester Matrix
    Sk = BuildSylvesterMatrix_3Polys(...
        lambda .* fww, ...
        mu .* gww, ...
        rho .* hww , m, n, o, k);
    
    % Get vector of singular values
    vSingularValues = svd(Sk);
    arr_SingularValues{i} = vSingularValues;
    arr_NormalisedSingularValues{i} = vSingularValues./vSingularValues(1);
    
    
    
    [~,R] = qr(Sk);
    [~,c] = size(R);
    
    arr_R1{i} = R(1:c, 1:c);
    
end



switch SETTINGS.RANK_REVEALING_METRIC
    % R1 Row Norms
    % R1 Row Diagonals
    % Singular Values
    % Residuals
    case 'R1 Row Norms'
        
        arr_R1_RowNorms = cell(nSubresultants,1);
        vMaxRowNormR1 = zeros(nSubresultants,1);
        vMinRowNormR1 = zeros(nSubresultants,1);
        
        for i = 1:1:nSubresultants
            
            arr_R1_RowNorms{i} = sqrt(sum(arr_R1{i}.^2,2))./norm(arr_R1{i});
            
            % Get maximum and minimum row norms of rows of R1.
            vMaxRowNormR1(i) = max(arr_R1_RowNorms{i});
            vMinRowNormR1(i) = min(arr_R1_RowNorms{i});
            
        end
        
        if(SETTINGS.PLOT_GRAPHS_RANK)
            
            plotRowNorms(arr_R1RowNorms, limits_k, limits_t)
            plotMaxMinRowNorms(vMaxRowNormR1, vMinRowNormR1, limits_k, limits_t, rank_range)
            
        end
        
        vMetric = log10(vMinRowNormR1 ./ vMaxRowNormR1);
        
    case 'R1 Row Diagonals'
        
        vMaxDiagonalEntry = zeros(nSubresultants, 1);
        vMinDiagonalEntry = zeros(nSubresultants, 1);
        
        for i = 1:1:length(arr_R1)
            
            % Get maximum diagonal
            vMaxDiagonalEntry(i) = max(abs(diag(arr_R1{i})));
            
            % Get minimum diagonal
            vMinDiagonalEntry(i) = min(abs(diag(arr_R1{i})));
            
        end
        
        vRatio_MaxMin_DiagonalEntry = vMinDiagonalEntry ./ vMaxDiagonalEntry;
        
        if(SETTINGS.PLOT_GRAPHS_RANK)
            
            %plotDiagonalsR1(arr_R1, myLimits, limits_t)
            plotMaxMinDiagonalR1(vRatio_MaxMin_DiagonalEntry, limits_k, limits_t, rank_range);
            
        end
        
        vMetric = log10(vRatio_MaxMin_DiagonalEntry);
        
    case 'Minimum Singular Values'
        
        vMinimumSingularValues = zeros(length(arr_SingularValues),1);
        
        for i = 1:1:length(arr_SingularValues)
            
            vSingularValues = arr_SingularValues{i};
            
            vMinimumSingularValues(i) = min(vSingularValues);
            
        end
        
        if(SETTINGS.PLOT_GRAPHS_RANK)
            plotSingularValues_S1(arr_SingularValues{1})
            plotSingularValues(arr_SingularValues, limits_k, limits_t);
            plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t, rank_range);
            
            
        end
        
        vMetric = log10(vMinimumSingularValues);
        
    case 'Normalised Minimum Singular Values'
        
        vMinimumNormalisedSingularValues = zeros(length(arr_NormalisedSingularValues),1);
        
        
        for i = 1:1:length(arr_NormalisedSingularValues)
            
            % Get set of all normalised singular values of S_{i}
            vNormalisedSingularValues = arr_NormalisedSingularValues{i};
            
            % Get minimum of the set
            vMinimumNormalisedSingularValues(i) = min(vNormalisedSingularValues);
            
        end
        
        if(SETTINGS.PLOT_GRAPHS_RANK)
            
            
            plotSingularValues(arr_NormalisedSingularValues, limits_k, limits_t);
            plotMinimumSingularValues(vMinimumNormalisedSingularValues, limits_k, limits_t, rank_range);
            
        end
        
        vMetric = log10(vMinimumNormalisedSingularValues);
        
    case 'Residuals'
        
        error('err')
        
    otherwise
        error('err')
end


if (lowerLimit_k == upperLimit_k)
    
    % Get the degree of the GCD
    t = GetGCDDegree_OneSubresultant(vMetric);
    
    GM_fx = vGM_fx(t - lowerLimit_k + 1);
    GM_gx = vGM_gx(t - lowerLimit_k + 1);
    GM_hx = vGM_hx(t - lowerLimit_k + 1);
    
    lambda = vAlpha(t - lowerLimit_k + 1);
    mu = vBeta(t - lowerLimit_k + 1);
    rho = vGamma(t - lowerLimit_k + 1);
    th1 = vTh1(t - lowerLimit_k + 1);
    th2 = vTh2(t - lowerLimit_k + 1);
    
else
    
    % Get the degree of the GCD
    t = GetGCDDegree_MultipleSubresultants(vMetric, limits_k, limits_t, rank_range);
    
    try
    rank_range_low = vMetric(t - (lowerLimit_k - 1) );
    if (t < upperLimit_k)
        rank_range_high = vMetric(t - (lowerLimit_k - 1) + 1);
    end
    catch
        fprintf('out of range')
    end
    
    % update rank range
    rank_range = [rank_range_low rank_range_high];
    
    
    GM_fx = vGM_fx(t - lowerLimit_k + 1);
    GM_gx = vGM_gx(t - lowerLimit_k + 1);
    GM_hx = vGM_hx(t - lowerLimit_k + 1);
    
    lambda = vAlpha(t - lowerLimit_k + 1);
    mu = vBeta(t - lowerLimit_k + 1);
    rho = vGamma(t - lowerLimit_k + 1);
    th1 = vTh1(t - lowerLimit_k + 1);
    th2 = vTh2(t - lowerLimit_k + 1);
    
end

end


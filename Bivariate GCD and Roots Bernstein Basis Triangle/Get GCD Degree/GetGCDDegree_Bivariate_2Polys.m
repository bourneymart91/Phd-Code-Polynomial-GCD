function [t, GM_fx, GM_gx, alpha, th1, th2, rank_range] = GetGCDDegree_Bivariate_2Polys(fxy, gxy, m, n, limits_t, rank_range)
% GetGCDDegree(fxy, gxy, m, n)
%
% Inputs.
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% limits_t : (Int Int)
%
% rank_range : [Float Float]
%
% % Outputs
%
% t : (Int) Degree of GCD d(x,y)
%
% GM_fx : (Float) Geometric Mean of f(x,y)
%
% GM_gx : (Float) Geometric Mean of g(x,y)
%
% alpha : (Float) Optimal value of \alpha
%
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}

global SETTINGS

% Set limits for k
limits_k = [1 min(m,n)];
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

% Get number of Sylvester subresultants
nSubresultants = upperLimit_k - lowerLimit_k + 1;

% Initialise some arrays
arr_SingularValues = cell(nSubresultants,1);
arr_NormalisedSingularValues = cell(nSubresultants,1);

arr_R1 = cell(nSubresultants,1);

vGM_fx = zeros(nSubresultants,1);
vGM_gx = zeros(nSubresultants,1);

vAlpha = zeros(nSubresultants,1);
vTh1 = zeros(nSubresultants,1);
vTh2 = zeros(nSubresultants,1);

for i = 1:1:nSubresultants
    
    k = lowerLimit_k + (i-1);
    
    % %
    % Geometric Mean
    
    % Get Geometric mean of f(x,y) in C_{n-k}(f) and g(x,y) in C_{m-k}(g)
    vGM_fx(i) = GetMean(fxy, m, n-k);
    vGM_gx(i) = GetMean(gxy, n, m-k);
    
    % Divide entries of f(x,y) and g(x,y) by geometric mean
    fxy_n = fxy ./ vGM_fx(i);
    gxy_n = gxy ./ vGM_gx(i);
    
    % %
    % Preprocess
    [vAlpha(i), vTh1(i), vTh2(i)] = Preprocess(fxy_n, gxy_n, m, n, k);
    
    % Get f(x,y) and g(x,y) with thetas to get f(w_{1},w_{2}) and g(w_{1},w_{2})
    fww = GetWithThetas(fxy_n, m, vTh1(i), vTh2(i));
    gww = GetWithThetas(gxy_n, n, vTh1(i), vTh2(i));
    
    
    % Plot coefficients of unprocessed and preprocessed polynomails (Only on
    % first iteration, to avoid too many graphs
    if i == 1
        if SETTINGS.PLOT_GRAPHS_PREPROCESSING
            
            PlotPolynomials(...
                {fxy, fww, gxy, vAlpha(i)*gww},...
                {...
                '$f\left(x, y \right)$',...
                '$\tilde{f}\left(\omega_{1},\omega_{2}\right)$',...
                'g(x,y)',...
                '$\alpha \tilde{g} \left(\omega_{1}, \omega_{2} \right)$'...
                });
        end
        
    end
    
    % %
    % Build the Sylvester Matrix
    Sk = BuildSylvesterMatrix_2Polys(fww, vAlpha(i).*gww, m, n, k);
    
    % Get vector of singular values
    
    vSingularValues = svd(Sk);
    
    arr_SingularValues{i} = vSingularValues;
    arr_NormalisedSingularValues{i} = vSingularValues./vSingularValues(1);
    
    
    
    
    [~,R] = qr(Sk);
    [~,c] = size(R);
    
    arr_R1{i} = R(1:c, 1:c);
    
end

% Metric used to compute the degree of the GCD
% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals
switch SETTINGS.RANK_REVEALING_METRIC
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
        
        if (SETTINGS.PLOT_GRAPHS_RANK)
            
            plotRowNorms(arr_R1_RowNorms, limits_k, limits_t)
            plotMaxMinRowNorms(vMaxRowNormR1, vMinRowNormR1, limits_k, limits_t, rank_range)
            
        end
        
        vMetric = log10(vMinRowNormR1 ./ vMaxRowNormR1);
        
    case 'R1 Row Diagonals'
        
        vMaxDiagonalEntry = zeros(nSubresultants,1);
        vMinDiagonalEntry = zeros(nSubresultants,1);
        
        for i = 1:1:length(arr_R1)
            
            % Get maximum diagonal
            vMaxDiagonalEntry(i) = max(abs(diag(arr_R1{i})));
            
            % Get minimum diagonal
            vMinDiagonalEntry(i) = min(abs(diag(arr_R1{i})));
            
        end
        
        vRatio_MaxMin_DiagonalEntry = vMinDiagonalEntry ./ vMaxDiagonalEntry;
        
        if(SETTINGS.PLOT_GRAPHS_RANK)
            %plotDiagonalsR1(arr_R1, limits_k, limits_t)
            plotMaxMinDiagonalR1(vRatio_MaxMin_DiagonalEntry, limits_k, limits_t, rank_range);
        end
        
        vMetric = log10(vRatio_MaxMin_DiagonalEntry);
        
    case 'Normalised Minimum Singular Values'
        
        
        vMinimumNormalisedSingularValues = zeros(length(arr_SingularValues),1);
        
        
        for i = 1:1:length(arr_SingularValues)
            
            vSingularValues = arr_SingularValues{i};
            vMinimumSingularValues(i) = min(vSingularValues);
            
            vSingularValues = arr_NormalisedSingularValues{i};
            vMinimumNormalisedSingularValues(i) = min(vSingularValues);
            
        end
        
        if(SETTINGS.PLOT_GRAPHS_RANK)
            
            
            plotSingularValues(arr_NormalisedSingularValues, limits_k, limits_t);
            %plotSingularValues(arr_SingularValues, limits_k, limits_t);
            
            %plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t, rank_range);
            plotMinimumSingularValues(vMinimumNormalisedSingularValues, limits_k, limits_t, rank_range);
        end
        
        %vMetric = log10(vMinimumSingularValues);
        vMetric = log10(vMinimumNormalisedSingularValues);
        
        
    case 'Minimum Singular Values'
        
        vMinimumSingularValues = zeros(length(arr_SingularValues),1);
        
        for i = 1:1:length(arr_SingularValues)
            
            vSingularValues = arr_SingularValues{i};
            vMinimumSingularValues(i) = min(vSingularValues);
            
            
        end
        
        if(SETTINGS.PLOT_GRAPHS_RANK)
            
            plotSingularValues(arr_SingularValues, limits_k, limits_t);
            
            plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t, rank_range);
            
        end
        
        vMetric = log10(vMinimumSingularValues);
        
    case 'Residuals'
        error('err')
        
    otherwise
        error('Error : Not a valid metric')
end

if (lowerLimit_k == upperLimit_k)
    
    t = GetGCDDegree_OneSubresultant(vMetric);
    
else
    
    t = GetGCDDegree_MultipleSubresultants(vMetric,[lowerLimit_k, upperLimit_k], limits_t, rank_range);
    
    fprintf('Computed Value t : %i',t);
    
    rank_range_low = vMetric(t - lowerLimit_k) + 1;
    rank_range_high = vMetric(t - lowerLimit_k) + 1;
    rank_range = [rank_range_low rank_range_high];
    
    
    
    LineBreakMedium()
    LineBreakMedium()
    
    GM_fx = vGM_fx(t-lowerLimit_k + 1);
    GM_gx = vGM_gx(t-lowerLimit_k + 1);
    alpha = vAlpha(t-lowerLimit_k + 1);
    th1 = vTh1(t-lowerLimit_k + 1);
    th2 = vTh2(t-lowerLimit_k + 1);
    
end

end




function PlotPolynomials(arrPolys, arrPolyNames)
%
% % Inputs
%
% fxy : (matrix)
%
% fww : (matrix) Coefficients of preprocessed polynomial


% Get number of polynomails
nPolys = length(arrPolys);


figure()

hold on

for i = 1 : 1 : nPolys
    
    
    fxy = arrPolys{i};
    
    % Get degree of polynomial
    m = GetDegree_Bivariate(fxy);
    
    % Get number of coefficients
    nCoefficients_fxy = nchoosek(m + 2, 2);
    
    % Get as vector
    v_fxy = GetAsVector(fxy);
    v_fxy = v_fxy(1 : nCoefficients_fxy);
    
    % Get polynomial name
    poly_name = arrPolyNames{i};
    
    % Plot coefficients
    plot(1 : 1 : nCoefficients_fxy, log10(v_fxy), '-o', 'DisplayName', poly_name);
    
    
end


% axis labels
%
xlabel('$i$ : Coefficient Index','Interpreter','latex','FontSize', 20)
ylabel('$\log_{10} \left( \Re \right)$','Interpreter', 'latex', 'FontSize', 20)

% Legend settings
leg1 = legend(gca,'show');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',20);
hold off


% Resizing Figure
myplot = gca;
myval_side = 0.12;
myval_base = 0.10;
set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])

% Set window size
set(gcf, 'Position', [100, 100, 700, 700])

% Display Properties
grid on
box on

end
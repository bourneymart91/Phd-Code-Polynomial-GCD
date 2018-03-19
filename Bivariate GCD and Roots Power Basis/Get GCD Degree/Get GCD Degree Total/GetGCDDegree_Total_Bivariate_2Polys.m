function [t] = GetGCDDegree_Total_Bivariate_2Polys(fxy, gxy, m, n, limits_t, rank_range)
% Calculate the degree of the GCD of two bivariate power basis polynomials.
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)  
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% m : (Int) Total degree of polynomial f(x,y) 
% 
% n : (Int) Total degree of polynomial g(x,y)
%
% limits_t : [Int Int] 
%
% rank_range : [Float Float]
%
%
% % Outputs
%
% t : (Int) Total degree of the GCD d(x,y)

% Set upper and lower limits
limits_k = [1, min(m,n)];
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

% Get number of Sylvester subresultant matrices
nSubresultants = upperLimit_k - lowerLimit_k + 1;

% Initialise some cell arrays
arr_R1_RowNorm = cell(nSubresultants, 1);
arr_R1_Diag = cell(nSubresultants, 1);
arr_SingularValues = cell(nSubresultants, 1);

% Pad the coefficients of f(x,y) and g(x,y0
% this is equivalent to degree elevating so that f is of degree (m,m), and
% g is of degree (n,n)
fxy_matrix_padd = zeros(m+1, m+1);
gxy_matrix_padd = zeros(n+1, n+1);


[m1, m2] = GetDegree_Bivariate(fxy);
fxy_matrix_padd(1:m1+1 , 1:m2+1) = fxy;

[n1, n2] = GetDegree_Bivariate(gxy);
gxy_matrix_padd(1:n1+1 , 1:n2+1) = gxy;



% let k represent the total degree of the common divisor d_{k}(x,y)
for i = 1 : 1 : nSubresultants
    
    k = lowerLimit_k + (i-1) ;
    

    % Build the partitions of the Sylvester matrix
    T_f = BuildT1_Total_Bivariate(fxy_matrix_padd, m, n-k);
    T_g = BuildT1_Total_Bivariate(gxy_matrix_padd, n, m-k);


    % Build the sylvester matrix
    Sk = [T_f T_g];
        
    
    % Using QR Decomposition of the sylvester matrix
    [~,R] = qr(Sk);
    
    % Take absolute values.
    R = abs(R);
    
    % Get number of rows in R1
    [R1_rows,~] = size(diag(R));
    
    % Obtain R1 the top square of the R matrix.
    R1 = R(1:R1_rows,1:R1_rows);
    
    % Get Norms of each row in the matrix R1
    arr_R1_RowNorm{i} = sqrt(sum(R1.^2,2))./norm(R1);
    
    % Get ONLY the diagonal elements and normalise them.
    arr_R1_Diag{i} = diag(R1);
    
    % Get singular values of S_{k}
    arr_SingularValues{i} = svd(Sk);
    
end



% R1 Row Norms
% R1 Row Diagonals
% Minimum Singular Values
% Residuals

global SETTINGS

switch SETTINGS.RANK_REVEALING_METRIC
    
    case 'Minimum Singular Values'
        
        % Initialise a vector to store minimum singular values
        vMinimumSingularValue = zeros(nSubresultants,1);
        
        for i = 1:1:nSubresultants
            vMinimumSingularValue(i) = min(arr_SingularValues{i});
        end
        
        if (SETTINGS.PLOT_GRAPHS)
            
            plotMinimumSingularValues_degreeTotal(vMinimumSingularValue, limits_k, limits_t, rank_range);
            plotSingularValues_degreeTotal(arr_SingularValues, limits_k, limits_t);
            
        end
        
        vMetric = log10(vMinimumSingularValue);
        
    case 'R1 Row Norms'
        % Get max/min row norms
        
        % Initialise vectors to store maximum row norm and minimum row norm
        vMaxRowNormR1 = zeros(nSubresultants,1);
        vMinRowNormR1 = zeros(nSubresultants,1);
        
        for i = 1:1:nSubresultants
            % Get the maximum row norm
            vMaxRowNormR1(i) = max(arr_R1_RowNorm{i});
            
            % Get the minimal row norm
            vMinRowNormR1(i) = min(arr_R1_RowNorm{i});
        end
        
        vRatio_MaxMin_RowNorm_R1 = vMinRowNormR1./vMaxRowNormR1;
        
        
        
        if (SETTINGS.PLOT_GRAPHS)
            
            plotRowNorm_degreeTotal(arr_R1_RowNorm, limits_k, limits_t, rank_range)
            plotMaxMinRowNorm_degreeTotal(vRatio_MaxMin_RowNorm_R1, limits_k, limits_t, rank_range);
            
        end
        
        vMetric = log10(vRatio_MaxMin_RowNorm_R1);
        
    case 'R1 Row Diagonals'
        
        % Initialise vectors
        v_maxDiagR1 = zeros(nSubresultants,1);
        v_minDiagR1 = zeros(nSubresultants,1);
        
        for i = 1:1:nSubresultants
            % Get the maximum diagonal of R1
            v_maxDiagR1(i) = max(abs(arr_R1_Diag{i}));
            
            % Get the minimum diagonal of R1
            v_minDiagR1(i) = min(abs(arr_R1_Diag{i}));
        end
        
        
        % Get max/min diagonal entries
        vRatio_MaxMin_Diags_R1 = v_minDiagR1./v_maxDiagR1;
        
        
        if(SETTINGS.PLOT_GRAPHS)
            plotRowDiag_degreeTotal(arr_R1_Diag, limits_k, limits_t);
            plotMaxMinRowDiag_degreeTotal(vRatio_MaxMin_Diags_R1, limits_k, limits_t);
        end
        
        vMetric = log10(vRatio_MaxMin_Diags_R1);
        
    case 'Residuals'
        
        error('Not Developed')
        
    otherwise
        error('%s is not a valid method',SETTINGS.RANK_REVEALING_METRIC);
        
end


if i == 1
    fprintf('Only one subresultant exists')
end


% if only one subresultant exists.
if lowerLimit_k == upperLimit_k
    t = GetGCDDegree_OneSubresultant(Sk);
    return;
else
    t = GetGCDDegree_MultipleSubresultants(vMetric, limits_k, limits_t, rank_range );
end

% Plot graphs
%PlotGraphs()

end

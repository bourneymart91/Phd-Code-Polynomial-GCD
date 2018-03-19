function [t, rank_range] = GetGCDDegree_Total_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, limits_t, rank_range)
% Calculate the degree of the GCD of two bivariate Power Basis polynomials.
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of the polynomial h(x,y)
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% o : (Int) Total degree of h(x,y)
%
% limits_t : (Int Int) Minimum and Maximum values bounding the computation of t
%
% rank_range : [Float Float]
%
% % Outputs
%
% t : (Int) Total degree of the GCD d(x,y).


limits_k = [0 min([m n o])];

% Set upper and lower bound for total degree t.
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

rank_range_low = rank_range(1);
rank_range_high = rank_range(2);

% if the upper and lower bound are equal, and not equal to one, then set
% the total degree to lower bound. If bound = 1, then it is possible for
% the polynomials to be coprime.
if (lowerLimit_k == upperLimit_k && lowerLimit_k ~= 1)
    
    t = lowerLimit_k;
    return;
    
end



% %
% pad the coefficients of fxy and gxy
% this is equivalent to degree elevating so that f is of degree (m,m), and
% g is of degree (n,n)
fxy_matrix_padd = zeros(m+1, m+1);
gxy_matrix_padd = zeros(n+1, n+1);
hxy_matrix_padd = zeros(o+1, o+1);

% Get the degree of f(x,y), g(x,y) and h(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);

%
fxy_matrix_padd(1:m1+1, 1:m2+1) = fxy;
gxy_matrix_padd(1:n1+1, 1:n2+1) = gxy;
hxy_matrix_padd(1:o1+1, 1:o2+1) = hxy;

% Get the number of subresultants
nSubresultants = upperLimit_k - lowerLimit_k + 1;

% Initialise some cell arrays
arr_R1_RowNorm = cell(nSubresultants, 1);
arr_R1_Diag = cell(nSubresultants, 1);
arr_SingularValues = cell(nSubresultants, 1);


% let k represent the total degree of the common divisor
for i = 1 : 1 : nSubresultants
    
    k = lowerLimit_k + (i-1);
    
    % Build the partitions T1 and T2 of the Sylvester matrix
    T1 = BuildT1_Total_Bivariate(fxy_matrix_padd, m, n-k);
    T2 = BuildT1_Total_Bivariate(fxy_matrix_padd, m, o-k);
    T3 = BuildT1_Total_Bivariate(gxy_matrix_padd, n, m-k);
    T4 = BuildT1_Total_Bivariate(hxy_matrix_padd, o, m-k);
    
    diagonal = blkdiag(T1,T2);
    column = [T3;T4];
    
    % Build the sylvester matrix
    Sk = [diagonal column];
    
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
            
            plotSingularValues_degreeTotal(arr_SingularValues, limits_k, limits_t);
            plotMinimumSingularValues_degreeTotal(vMinimumSingularValue, limits_k, limits_t, rank_range);
            
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
        
        vMetric = log10(vRatio_MaxMin_RowNorm_R1);
        
        if (SETTINGS.PLOT_GRAPHS)
            plotRowNorm_degreeTotal(arr_R1_RowNorm, limits_k, limits_t)
            plotMaxMinRowNorm_degreeTotal(vRatio_MaxMin_RowNorm_R1, limits_k, limits_t);
        end
        
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
        vMetric = log10(vRatio_MaxMin_Diags_R1);
        
        if(SETTINGS.PLOT_GRAPHS)
            plotRowDiag_degreeTotal(arr_R1_Diag, limits_k, limits_t);
            plotMaxMinRowDiag_degreeTotal(vRatio_MaxMin_Diags_R1, limits_k, limits_t);
        end
        
    case 'Residuals'
        error('Not Developed')
        
    otherwise
        
        error('%s is not a valid rank revealing metric', SETTINGS.RANK_REVEALING_METRIC)
end


if i == 1
    fprintf('Only one subresultant exists')
end


% if only one subresultant exists.
if lowerLimit_k == upperLimit_k
    
    t = GetGCDDegree_OneSubresultant(Sk);
    return;
else if lowerLimit_t == upperLimit_t
        t = lowerLimit_t;
    return;
else
    
    t = GetGCDDegree_MultipleSubresultants(vMetric, limits_k, limits_t, rank_range);
    
end

i = t - lowerLimit_k + 1;

% Define new rank range
rank_range_low = vMetric(i);
if i == upperLimit_k - lowerLimit_k + 1
else
    rank_range_high = vMetric(i+1);
end

rank_range = [rank_range_low rank_range_high];


end



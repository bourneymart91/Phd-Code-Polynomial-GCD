function [t1,t2] = GetGCDDegree_Relative_Bivariate_2Polys_Fast(fxy, gxy, limits_t1, limits_t2)
% GetGCDDegree_Relative_2Polys(fxy,gxy)
%
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two
% polynomials f(x,y) and g(x,y)
%
% Inputs.
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% limits_t1 : (Int Int) 
%
% limits_t2 : (Int Int)
%
% % Outputs
%
% t1 : (Int) Degree of d(x,y) with respect to x
%
% t2 : (Int) Degree of d(x,y) with respect to y


% Get degree of f(x,y) and g(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

% Get my limits
limits_k1 = [0 min(m1, n1)];
limits_k2 = [0 min(m2, n2)];

lowerLimit_t1 = limits_k1(1);
upperLimit_t1 = limits_k1(2);
lowerLimit_t2 = limits_k2(1);
upperLimit_t2 = limits_k2(2);

% Get number of Sylvester subresultants to be constructed
nSubresultants_k1 = upperLimit_t1 - lowerLimit_t1 + 1;
nSubresultants_k2 = upperLimit_t2 - lowerLimit_t2 + 1;

% Initialise cell array to store all singular values
arr_SingularValues = cell( nSubresultants_k1, nSubresultants_k2 );

% Initialise cell array to store all R matrices from QR decomposition of
% each S_{k_{1},k_{2}}
arr_R1 = cell(nSubresultants_k1, nSubresultants_k2);
arr_R1_RowNorms = cell(nSubresultants_k1, nSubresultants_k2);


% Algorithm for producing set of S_{k1,k2} sylvester matrices
% Produce S_{i,0} as polynomials in terms of x
% Produce S_{0,j} for j = 1,...,min(m2,n2) by removing rows and column
% sections where Sylvester matrices are built in terms of polynomials in x.

arr_Sk1k2 = cell(nSubresultants_k1,nSubresultants_k2);
arr_Q = cell(nSubresultants_k1,nSubresultants_k2);
arr_R = cell(nSubresultants_k1,nSubresultants_k2);

% For each Sylvester subresultant
for i1 = 1:1:nSubresultants_k1
    
    k1 = lowerLimit_t1 + (i1-1);
    
    % Produce first row of Sylvester matrices.
    arr_Sk1k2{i1,1} = BuildT_Relative_Bivariate_2Polys(fxy, gxy, k1, 0);
    
    % Get QR decomposition
    [arr_Q{i1,1}, arr_R{i1,1}] = qr(arr_Sk1k2{i1,1});
    [~,c] = size(arr_R{i1,1});
    
    % Get matrix R_{1}
    arr_R1{i1,1} = arr_R{i1,1}(1:c,1:c);
    
    
    arr_SingularValues{i1, 1} = svd(arr_Sk1k2{i1,1});
    
    for i2 = 2:1:nSubresultants_k2
        
        k2 = lowerLimit_t2 + (i2-1);
        
        % Get the number of column sections in the previous
        nColumnPartitions_T1_prev = (n2-(k2-1)+1);
        nColumnPartitions_T1_new = (n2-k2+1);
        
        % Each section has
        nColumnsPerSection_T1 = (n1-k1+1);
        
        nColumns_T1_prev = nColumnPartitions_T1_prev * nColumnsPerSection_T1;
        nColumns_T1_new = nColumnPartitions_T1_new * nColumnsPerSection_T1;
        
        nColsToBeRemoved_T1 = nColumns_T1_prev - nColumns_T1_new;
        
        % %
        %
        
        nColumnPartitions_T2_prev = (m2-(k2-1)+1);
        nColumnPartitions_T2_new = (m2-k2+1);
        
        nColumnsPerSection_T2 = (m1-k1+1);
        
        nColumns_T2_prev = nColumnPartitions_T2_prev * nColumnsPerSection_T2;
        nColumns_T2_new = nColumnPartitions_T2_new * nColumnsPerSection_T2;
        
        nColsToBeRemoved_T2 = nColumns_T2_prev - nColumns_T2_new;
        
        nColumns_Sk1k2_prev = nColumns_T1_prev + nColumns_T2_prev;
                
        nColsToBeRemoved = nColsToBeRemoved_T1 + nColsToBeRemoved_T2;
        
        % %
        %
        nRowPartitions_prev = m2 + n2 - (k2-1) + 1;
        nRowPartitions_new = m2 + n2 - k2 + 1;
        
        nRowsPerPartition = m1+n1-k1+1;
        
        nRows_prev = nRowPartitions_prev * nRowsPerPartition;
        nRows_new = nRowPartitions_new * nRowsPerPartition;
        
        nRowsToBeRemoved = nRows_prev - nRows_new;
        
        % %
        % Get index of columns to be removed
        idx_lastColDelete_T1 = nColumns_T1_prev;
        idx_firstColDelete_T1 = idx_lastColDelete_T1 - nColsToBeRemoved_T1 + 1;
        
        idx_vec_T1 = idx_firstColDelete_T1 : 1 : idx_lastColDelete_T1;
        
        idx_lastColDelete_T2 = nColumns_Sk1k2_prev;
        idx_firstColDelete_T2 = idx_lastColDelete_T2 - nColsToBeRemoved_T2 + 1;
        
        idx_vec_T2 = idx_firstColDelete_T2 : 1 : idx_lastColDelete_T2;
        
        vec_cols = [idx_vec_T1 idx_vec_T2];
        
        % Get index of rows to be removed
        
        idx_lastRowDelete = nRows_prev;
        idx_firstRowDelete = nRows_prev - nRowsToBeRemoved + 1;
        vec_rows = idx_firstRowDelete : 1 : idx_lastRowDelete;
        
        %%
        
        Sk1k2_prev = arr_Sk1k2{i1,i2-1};
        
        % Remove rows from previous sylvester matrix
        Sk1k2_prev(vec_rows,:) = [];
        
        % Remove columns from previous sylvester matrix
        Sk1k2_prev(:,vec_cols) = [];
        
        Sk1k2_test = BuildT_Relative_Bivariate_2Polys(fxy, gxy, k1, k2);
        
        arr_Sk1k2{i1,i2} = Sk1k2_prev;
        
        % Get the singular values of S_{k_{1},k_{2}}
        arr_SingularValues{i1, i2} = svd(arr_Sk1k2{i1,i2});
        
        %[~, R_temp] = qr(Sk1k2_test);
        
        %%
        % Remove all rows and columns
        
        Q_temp = arr_Q{i1,i2-1}; 
        R_temp = arr_R{i1,i2-1};
        
        for i = 1:1:nColsToBeRemoved

            temp_idx = vec_cols(end);

            vec_cols(end) = [];

            [Q_temp, R_temp] = qrdelete(Q_temp, R_temp, temp_idx,'col');

        end
        
        
        
        for i = 1:1:nRowsToBeRemoved
            [Q_temp, R_temp] = qrdelete(Q_temp,R_temp, idx_lastRowDelete - (i-1),'row');
        end
        
                
        arr_Q{i1,i2} = Q_temp;
        arr_R{i1,i2} = R_temp;
        
        [~,c] = size(arr_R{i1,i2});
        R1 = arr_R{i1,i2}(1:c,1:c);
        arr_R1{i1,i2} = R1;
        
        
       
        
        
        
    end
end


global SETTINGS
% Metric used to compute the degree of the GCD
% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals
fprintf('Metric used to compute degree of GCD : %s \n', SETTINGS.RANK_REVEALING_METRIC)
switch SETTINGS.RANK_REVEALING_METRIC
    case 'Minimum Singular Values'
        
        
        mat_MinimumSingularValues = zeros(nSubresultants_k1, nSubresultants_k2);
        for i1 = 1:1:nSubresultants_k1
            for i2 = 1:1:nSubresultants_k2
                
                %k1 = lowerLimit_k1 + (i1-1);
                %k2 = lowerLimit_k2 + (i2-1);
                
                mat_MinimumSingularValues(i1, i2) = min(arr_SingularValues{i1,i2});
            end
        end
        
        mat_metric = mat_MinimumSingularValues;
        
        if (SETTINGS.PLOT_GRAPHS)
            plotSingularValues_degreeRelative(arr_SingularValues, limits_k1, limits_k2, limits_t1, limits_t2)
            plotMinimumSingularValues_degreeRelative(mat_MinimumSingularValues, limits_k1, limits_k2, limits_t1, limits_t2);
        end
        
    case 'R1 Row Norms'
        matMaxRowNorm = zeros(nSubresultants_k1 , nSubresultants_k2);
        matMinRowNorm = zeros(nSubresultants_k1 , nSubresultants_k2);
        for i1 = 1:1:nSubresultants_k1
            for i2 = 1:1:nSubresultants_k2
                
                arr_R1_RowNorms{i1,i2} = sqrt(sum(arr_R1{i1, i2}.^2,2))./norm(arr_R1{i1, i2});
                
                matMaxRowNorm(i1,i2) = max(arr_R1_RowNorms{i1,i2});
                matMinRowNorm(i1,i2) = min(arr_R1_RowNorms{i1,i2});
                
            end
        end
        
        if (SETTINGS.PLOT_GRAPHS)
            
            plotRowNormsR1_degreeRelative(arr_R1_RowNorms, limits_k1, limits_k2, limits_t1, limits_t2);
            plotMaxMinRowNormsR1_degreeRelative(matMaxRowNorm, matMinRowNorm, limits_k1, limits_k2, limits_t1, limits_t2);
            
        end
        mat_metric = matMaxRowNorm ./ matMinRowNorm;
        
    case 'R1 Row Diagonals'
        
        mat_MinR1Diagonal = zeros(nSubresultants_k1, nSubresultants_k2);
        mat_MaxR1Diagonal = zeros(nSubresultants_k1, nSubresultants_k2);
        arr_DiagR1 = cell(nSubresultants_k1, nSubresultants_k2);
        
        for i1 = 1:1:nSubresultants_k1
            for i2 = 1:1:nSubresultants_k2
                
                
                
                arr_DiagR1{i1, i2} = diag(arr_R1{i1, i2});
                
                % Get matrix of maximum values on diagonals of R_{i1,i2}
                mat_MaxR1Diagonal(i1, i2) = max(abs(arr_DiagR1{i1, i2}));
                
                % Get matrix of minimum values on diagonals of R_{i1,i2}
                mat_MinR1Diagonal(i1, i2) = min(abs(arr_DiagR1{i1, i2}));
            end
        end
        
        mat_Ratio = mat_MinR1Diagonal ./ mat_MaxR1Diagonal;
        
        mat_metric = mat_Ratio;
        
        if (SETTINGS.PLOT_GRAPHS)
            plotMaxMinDiagonals_degreeRelative(mat_MaxR1Diagonal, mat_MinR1Diagonal, limits_k1, limits_k2, limits_t1, limits_t2);
            plotDiagonalsR1_degreeRelative(arr_DiagR1, limits_k1, limits_k2, limits_t1, limits_t2);
        end
    case 'Residuals'
        
        error('Code not yet developed for this branch')
        
    otherwise
        error('%s : is not a valid metric', SETTINGS.RANK_REVEALING_METRIC) 
end

% Compute the degree of the GCD
delta_x = diff(log10(mat_metric),1,1);
vec_delta_x = sum(delta_x,2);
[~, idx] = max(vec_delta_x);
t1 = lowerLimit_t1 + idx - 1;

delta_y = diff(log10(mat_metric),1,2);
vec_delta_y = sum(delta_y,1);
[~, idx] = max(vec_delta_y);
t2 = lowerLimit_t2 + idx - 1;



fprintf([mfilename ' : ' 'The Calculated Degree of the GCD is given by \n'])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt x : t1 = %i\n',t1)])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt y : t2 = %i\n',t2)])
end

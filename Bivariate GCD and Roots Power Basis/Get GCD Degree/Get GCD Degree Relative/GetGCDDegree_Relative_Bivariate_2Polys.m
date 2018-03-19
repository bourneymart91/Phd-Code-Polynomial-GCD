function [t1,t2] = GetGCDDegree_Relative_Bivariate_2Polys(fxy, gxy, limits_t1, limits_t2)
% GetGCDDegree_Relative_2Polys(fxy,gxy)
%
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two
% polynomials f(x,y) and g(x,y)
%
% Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
% 
% limits_t1 : (Int) (Int) : Lower and upper limit of possible value of t1
%
% limits_t2 : (Int) (Int) : Lower and upper limit of possible value of t2
%
% Outputs
%
% t1 : (Int) Degree of d(x,y) with respect to x
%
% t2 : (Int) Degree of d(x,y) with respect to y


% Get degree of f(x,y) and g(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

% Get limits
lowerLimit_t1 = limits_t1(1);
upperLimit_t1 = limits_t1(2);
lowerLimit_t2 = limits_t2(1);
upperLimit_t2 = limits_t2(2);

% Get limits in k1 and k_{2}
limits_k1 = [0 min(m1, n1)];
limits_k2 = [0 min(m2, n2)];

lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);
lowerLimit_k2 = limits_k2(1);
upperLimit_k2 = limits_k2(2);


% Get number of Sylvester subresultant matrices
nSubresultants_k1 = upperLimit_k1 - lowerLimit_k1 + 1;
nSubresultants_k2 = upperLimit_k2 - lowerLimit_k2 + 1;

% Initialise cell array to store all singular values
arr_SingularValues = cell( nSubresultants_k1, nSubresultants_k2 );

% Initialise cell array to store all R matrices from QR decomposition of
% each S_{k_{1},k_{2}}
arr_Sk1k2 = cell(nSubresultants_k1, nSubresultants_k2);
arr_R1 = cell(nSubresultants_k1, nSubresultants_k2);
arr_R1_RowNorms = cell(nSubresultants_k1, nSubresultants_k2);

% For each Sylvester subresultant
for i1 = 1:1:nSubresultants_k1
    for i2 = 1:1:nSubresultants_k2
        
        k1 = lowerLimit_t1 + (i1-1);
        k2 = lowerLimit_t2 + (i2-1);
        
        % Build the kth sylvester matrix
        Sk1k2 = BuildT_Relative_Bivariate_2Polys(fxy, gxy, k1, k2);
        arr_Sk1k2{i1,i2} = Sk1k2;
        
        % This test function is a rearranged for of the function above.
        % This form sees the polynomials f(x,y) and g(x,y) represented as a
        % by a set of polynomials only in x.
        %Sk1k2_test = BuildT_Relative_Bivariate_2Polys_NewMethod(fxy, gxy, k1, k2);
        %Sk1k2_test2 = BuildT_Relative_Bivariate_2Polys_NewMethod_alt(fxy, gxy, k1, k2);
        
        
        % Get the singular values of S_{k_{1},k_{2}}
        arr_SingularValues{i1, i2} = svd(Sk1k2);
        
        % Get R1 matrix
        %[~,R] = qr(Sk1k2);
        % Efficient QR when nRows > nCols. Does not explicitly compute Q
        % matrix.
        [R] = qr(Sk1k2,0);
        
        
        [~,c] = size(R);
        R1 = R(1:c,1:c);
        arr_R1{i1, i2} = R1;
        
        
        
    end
end


global SETTINGS
% Metric used to compute the degree of the GCD
% R1 Row Norms
% R1 Row Diagonals
% Minimum Singular Values
% Residuals

fprintf('Metric used to compute degree of GCD : %s', SETTINGS.RANK_REVEALING_METRIC);
switch SETTINGS.RANK_REVEALING_METRIC
    
    case 'Minimum Singular Values'
        
        % Initalise a matrix to store minimum Singular Values
        mat_MinimumSingularValues = zeros(nSubresultants_k1, nSubresultants_k2);
        
        for i1 = 1:1:nSubresultants_k1
            for i2 = 1:1:nSubresultants_k2
                               
                mat_MinimumSingularValues(i1,i2) = min(arr_SingularValues{i1,i2});
                
            end
        end
        
        mat_metric = log10(mat_MinimumSingularValues);
        
        % Plot Graphs
        if (SETTINGS.PLOT_GRAPHS)
            
            plotSingularValues_degreeRelative(arr_SingularValues, limits_k1, limits_k2, limits_t1, limits_t2)
            plotMinimumSingularValues_degreeRelative(mat_MinimumSingularValues, limits_k1, limits_k2, limits_t1, limits_t2);
            
        end
        
    case 'R1 Row Norms'
        
        % Initialise matrices to store max and minimum
        matMaxRowNorm = zeros(nSubresultants_k1 , nSubresultants_k2);
        matMinRowNorm = zeros(nSubresultants_k1 , nSubresultants_k2);
        
        for i1 = 1:1:nSubresultants_k1
            for i2 = 1:1:nSubresultants_k2
                
                arr_R1_RowNorms{i1, i2} = sqrt(sum(arr_R1{i1, i2}.^2,2))./norm(arr_R1{i1, i2});
                
                matMaxRowNorm(i1,i2) = max(arr_R1_RowNorms{i1,i2});
                matMinRowNorm(i1,i2) = min(arr_R1_RowNorms{i1,i2});
                
            end
        end
        
        % Plot Graphs
        if (SETTINGS.PLOT_GRAPHS)
            
            plotRowNormsR1_degreeRelative(arr_R1_RowNorms, limits_k1, limits_k2, limits_t1, limits_t2);
            plotMaxMinRowNormsR1_degreeRelative(matMaxRowNorm, matMinRowNorm, limits_k1, limits_k2, limits_t1, limits_t2);
            
        end
        
        mat_metric = log10(matMaxRowNorm ./ matMinRowNorm);
        
    case 'R1 Row Diagonals'
        
        mat_MinR1Diagonal = zeros(nSubresultants_k1, nSubresultants_k2);
        mat_MaxR1Diagonal = zeros(nSubresultants_k1, nSubresultants_k2);
        arr_DiagR1 = cell(nSubresultants_k1, nSubresultants_k2);
        
        for i1 = 1:1:nSubresultants_k1
            for i2 = 1:1:nSubresultants_k2
                
                
                
                arr_DiagR1{i1, i2} = diag(abs(arr_R1{i1, i2}));
                
                % Get matrix of maximum values on diagonals of R_{i1,i2}
                mat_MaxR1Diagonal(i1, i2) = max(arr_DiagR1{i1, i2});
                
                % Get matrix of minimum values on diagonals of R_{i1,i2}
                mat_MinR1Diagonal(i1, i2) = min(arr_DiagR1{i1, i2});
            end
        end
        
        mat_Ratio = mat_MinR1Diagonal ./ mat_MaxR1Diagonal;
        mat_metric = log10(mat_Ratio);
        
        % Plot Graphs
        if (SETTINGS.PLOT_GRAPHS)
            
            plotMaxMinDiagonals_degreeRelative(mat_MaxR1Diagonal, mat_MinR1Diagonal, limits_k1, limits_k2, limits_t1, limits_t2);
            plotDiagonalsR1_degreeRelative(arr_DiagR1, limits_k1, limits_k2, limits_t1, limits_t2);
            
        end
        
    case 'Residuals'
        
        error('Code not yet developed for this branch')
        
    otherwise
        error('%s : Not a valid metric', SETTINGS.RANK_REVEALING_METRIC)
end

% Compute the degree of the GCD
delta_x = diff((mat_metric),1,1);
vec_delta_x = sum(delta_x,2);
[~, idx] = max(vec_delta_x);
t1 = lowerLimit_t1 + idx - 1;

delta_y = diff((mat_metric),1,2);
vec_delta_y = sum(delta_y,1);
[~, idx] = max(vec_delta_y);
t2 = lowerLimit_t2 + idx - 1;



fprintf([mfilename ' : ' 'The Calculated Degree of the GCD is given by \n'])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt x : t1 = %i\n',t1)])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt y : t2 = %i\n',t2)])
end

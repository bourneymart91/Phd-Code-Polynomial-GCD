function [t1,t2] = GetGCDDegree_Relative_Bivariate_3Polys(fxy, gxy, hxy, limits_k1, limits_k2, limits_t1, limits_t2)
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two
% polynomials f(x,y) and g(x,y)
%
% Inputs.
%
% fxy : (Matrix) Coefficient matrix of polynomial f(x,y)
%
% gxy : (Matrix) Coefficient matrix of polynomial g(x,y)
%
% hxy : (Matrix) Coefficient matrix of polynomial h(x,y)
%
% myLimits_t1 : (Int) (Int)
%
% myLimits_t2 : (Int) (Int)
%
% limits_t1 : (Int) (Int)
% 
% limits_t2 : (Int) (Int)
%
% Outputs
%
% t1 : (Int) Degree of d(x,y) with respect to x
% 
% t2 : (Int) degree of d(x,y) with respect to y

% Get the degree structure of polynomial f(x,y), g(x,y) and h(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);

% Get limits
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);
lowerLimit_k2 = limits_k2(1);
upperLimit_k2 = limits_k2(2);

nSubresultants_k1 = upperLimit_k1 - lowerLimit_k1 + 1;
nSubresultants_k2 = upperLimit_k2 - lowerLimit_k2 + 1;


arr_SingularValues = cell(nSubresultants_k1, nSubresultants_k2);

% For each of the pairs [k1,k2]
for i1 = 1:1:nSubresultants_k1
    for i2 = 1:1:nSubresultants_k2
    
    k1 = lowerLimit_k1 + (i1-1);
    k2 = lowerLimit_k2 + (i2-1);
    
    % Build the partitions of the Sylvester matrix
    T1 = BuildT1_Relative_Bivariate(fxy, n1-k1, n2-k2);
    T2 = BuildT1_Relative_Bivariate(fxy, o1-k1, o2-k2);
    
    T3 = BuildT1_Relative_Bivariate(gxy, m1-k1, m2-k2);
    T4 = BuildT1_Relative_Bivariate(hxy, m1-k1, m2-k2);
    
    diagonal = blkdiag(T1,T2);
    column = [T3; T4];
    % Build the sylvester matrix
    Sk1k2 = [diagonal column];
    

   % Get the minimal singular value of matrix S_{k1,k2}
    arr_SingularValues{i1,i2} = svd(Sk1k2);
   
    
    
    
    end
end



global SETTINGS
% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals

switch SETTINGS.RANK_REVEALING_METRIC
    case 'Singular Values'
        
        mat_MinimumSingularValues = zeros( nSubresultants_k1 , nSubresultants_k2);
        
        for i1 = 1:1: nSubresultants_k1
            for i2 = 1:1: nSubresultants_k2
                
                %k1 = lowerLimit_k1 + (i1 - 1);
                %k2 = lowerLimit_k2 + (i2 - 1);
                
                mat_MinimumSingularValues(i1,i2) = min(arr_SingularValues{i1,i2});
                
            end
        end
        
        mat_metric = mat_MinimumSingularValues;
        
        %plotSingularValues_degreeRelative(arr_SingularValues, myLimits_t1, myLimits_t2, limits_t1, limits_t2);
        plotMinimumSingularValues_degreeRelative(mat_MinimumSingularValues, limits_k1, limits_k2, limits_t1, limits_t2);
        
    case 'R1 Row Norms'
        
        error('err')
        
    case 'R1 Row Diagonals'
        
        max_matrix = zeros( nSubresultants_k1, nSubresultants_k2);
        min_matrix = zeros( nSubresultants_k1, nSubresultants_k2);
        
        for i1 = 1:1: nSubresultants_k1
            for i2 = 1:1: nSubresultants_k2
                
                %k1 = lowerLimit_k1 + (i1-1);
                %k2 = lowerLimit_k2 + (i2-1);
                
                max_matrix(i1,i2) = max(abs(arr_DiagonalsR1{i1, i2}));
                min_matrix(i1,i2) = min(abs(arr_DiagonalsR1{i1, i2}));
                
                
            end
        end
        
        
        plotDiagonalsR1_degreeRelative(arr_DiagonalsR1, limits_k1, limits_k2, limits_t1, limits_t2);
        plotMaxMinDiagonals_degreeRelative(max_matrix,min_matrix, limits_k1, limits_k2, limits_t1, limits_t2);
        
        mat_metric = min_matrix./max_matrix;
        
        
        
    case 'Residuals'
        
        error('err');
        
    otherwise
        error('err');
end

% Compute the degree of the GCD
delta_x = diff(log10(mat_metric),1,1);
vec_delta_x = sum(delta_x,2);
[~, idx] = max(vec_delta_x);
t1 = lowerLimit_k1 + idx - 1;

delta_y = diff(log10(mat_metric),1,2);
vec_delta_y = sum(delta_y,1);
[~, idx] = max(vec_delta_y);
t2 = lowerLimit_k2 + idx - 1;


fprintf([mfilename ' : ' 'The Calculated Degree of the GCD is given by \n'])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt x : t1 = %i\n',t1)])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt y : t2 = %i\n',t2)])
end

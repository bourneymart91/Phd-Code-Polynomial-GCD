function [t1,t2] = GetGCDDegree_Relative_Given_t_Bivariate_2Polys(fxy, gxy, m, n, t, limits_t1, limits_t2)
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two
% polynomials f(x,y) and g(x,y)
%
% % Inputs.
%
% fxy : (Matrix) Coefficient matrix of polynomial f(x,y)
%
% gxy : (Matrix) Coefficient matrix of polynomial g(x,y)
%
% m : (Int) Total degree of polynomial f(x,y)
%
% n : (Int) Total degree of polynomial g(x,y)
%
% t : (Int) Total degree of GCD d(x,y)
%
% limits_t1 : [lowerLimit upperLimit]
%
% limits_t2 : [lowerLimit upperLimit]
%
% % Outputs
%
% t1 : (Int) Degree of d(x,y) with respect to x
%
% t2 : (Int) Degree of d(x,y) with respect to y

% Note this file differs from 'GetGCDDegreeRelative' since we make use of
% the computed value of 't' given t, for any k1,k2 pair, we know the
% location of certain zeros (if they exist in the lower right triangle) for
% polynomials u(x,y) and v(x,y). We remove the corresponding columns of
% S(f,g)

% Define my limits for the computation of the degree of the GCD. This can
% be set to include all subresultant matrices, or a subset defined by the
% computed limits passed in to this function.


[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

limits_k1 = [0 min(m1, n1)];
limits_k2 = [0 min(m2, n2)];


% Set upper and lower limits based on my limits.
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);
lowerLimit_k2 = limits_k2(1);
upperLimit_k2 = limits_k2(2);

% Get the number of Sylvester subresultant matrices to be computed
nSubresultants_k1 = upperLimit_k1 - lowerLimit_k1 + 1;
nSubresultants_k2 = upperLimit_k2 - lowerLimit_k2 + 1;

% Initialise arrays
arr_SingularValues = cell( nSubresultants_k1 , nSubresultants_k2);
arr_R = cell( nSubresultants_k1 , nSubresultants_k2);
arr_R1 = cell( nSubresultants_k1 , nSubresultants_k2);
arr_DiagonalsR1 = cell( nSubresultants_k1 , nSubresultants_k2);



% For each of the pairs [k1,k2]
for i1 = 1 : 1: nSubresultants_k1
    for i2 = 1 : 1: nSubresultants_k2
        
        k1 = lowerLimit_k1 + (i1 - 1);
        k2 = lowerLimit_k2 + (i2 - 1);
        
        % Build the Sylvester matrix S_{k,k1,k2}
        Skk1k2 = BuildT_Both_Bivariate_2Polys(fxy, gxy, m, n, t, k1, k2);
        
        
        % Get array of R matrices from QR decomposition of S_{k1,k2}
        [~,arr_R{i1,i2}] = qr(Skk1k2);
        [~,c] = size(arr_R{i1,i2});
        arr_R1{i1,i2} = arr_R{i1,i2}(1:c,1:c);
        arr_DiagonalsR1{i1,i2} = diag(arr_R1{i1,i2});
        
        
        % Get the minimal singular value of matrix S_{k1,k2}
        arr_SingularValues{i1,i2} = svd(Skk1k2);
        
        
    end
end




global SETTINGS
%   * R1 Row Norms
%   * R1 Row Diagonals
%   * Singular Values
%   * Residuals

switch SETTINGS.RANK_REVEALING_METRIC
    case 'Minimum Singular Values'
        
        mat_MinimumSingularValues = zeros( nSubresultants_k1 , nSubresultants_k2);
        
        for i1 = 1:1: nSubresultants_k1
            for i2 = 1:1: nSubresultants_k2
                
                %k1 = lowerLimit_k1 + (i1 - 1);
                %k2 = lowerLimit_k2 + (i2 - 1);
                
                mat_MinimumSingularValues(i1,i2) = min(arr_SingularValues{i1,i2});
                
            end
        end
        
        mat_metric = mat_MinimumSingularValues;
        
        plotSingularValues_degreeRelative(arr_SingularValues, limits_k1, limits_k2, limits_t1, limits_t2);
        plotMinimumSingularValues_degreeRelative(mat_MinimumSingularValues, limits_k1, limits_k2, limits_t1, limits_t2);
        
    case 'R1 Row Norms'
        
        error('err')
        
    case 'R1 Row Diagonals'
        
        max_matrix = zeros( nSubresultants_k1, nSubresultants_k2);
        min_matrix = zeros( nSubresultants_k1, nSubresultants_k2);
        
        for i1 = 1:1: nSubresultants_k1
            for i2 = 1:1: nSubresultants_k2
                

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

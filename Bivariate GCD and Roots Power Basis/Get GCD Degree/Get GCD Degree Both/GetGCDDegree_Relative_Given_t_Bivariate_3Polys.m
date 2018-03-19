function [t1,t2] = GetGCDDegree_Relative_Given_t_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, t, limits_t1, limits_t2, rank_range)
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two
% polynomials f(x,y) and g(x,y)
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of the polynomial h(x,y)
%
% m n o : (Int) (Int) (Int) Total degree of f(x,y), g(x,y) and h(x,y)
%
% t : (Int) Total degree of GCD d(x,y)
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

% Note this file differs from 'GetGCDDegreeRelative' since we make use of
% the computed value of 't' given t, for any k1,k2 pair, we know the
% location of certain zeros (if they exist in the lower right triangle) for
% polynomials u(x,y) and v(x,y). We remove the corresponding columns of
% S(f,g)

% Get the degree structure of polynomials f(x,y), g(x,y) and h(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);

% Set limits of k_{1} and k_{2}
limits_k1 = [0 min([m1, n1, o1])];
limits_k2 = [0 min([m2, n2, o2])];

% Get upper and lower limits of k_{1} and k_{2}
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);
lowerLimit_k2 = limits_k2(1);
upperLimit_k2 = limits_k2(2);

% Get number of Sylvester subresultants
nSubresultants_k1 = upperLimit_k1 - lowerLimit_k1 + 1;
nSubresultants_k2 = upperLimit_k2 - lowerLimit_k2 + 1;

arr_SingularValues = cell( nSubresultants_k1, nSubresultants_k2);
arr_R = cell(nSubresultants_k1, nSubresultants_k2);
arr_R1 = cell( nSubresultants_k1, nSubresultants_k2);
arr_Diagonals_R1 = cell( nSubresultants_k1, nSubresultants_k2);


% For each of the pairs [k1,k2]
for i1 = 1:1: nSubresultants_k1
    for i2 = 1:1: nSubresultants_k2
        
        
        k1 = lowerLimit_k1 + (i1 - 1);
        k2 = lowerLimit_k2 + (i2 - 1);
        
        
        % Build the Sylvester matrix S_{k,k1,k2}
        Skk1k2 = BuildT_Both_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, t, k1, k2);

        % Get the singular values of S_{k_{1},k_{2}}
        arr_SingularValues{i1,i2} = svd(Skk1k2);

    end
end



global SETTINGS
% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals

switch SETTINGS.RANK_REVEALING_METRIC
    case 'Minimum Singular Values'
        mat_MinimumSingularValues = zeros( nSubresultants_k1, nSubresultants_k2);
        
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
        
        error('Not Completed')
        
    case 'R1 Row Diagonals'
        
        error('Not Completed')
        
    case 'Residuals'
        
        error('Not Completed')
        
    otherwise 
        error(' %s is not a valid rank revealing metric', SETTINGS.RANK_REVEALING_METRIC)
end


% Compute the degree of the GCD
delta_x = diff(log10(mat_metric),1,1);
vDelta_x = sum(delta_x,2);
[maxDelta_x, idx] = max(vDelta_x);
% Determine whether the delta is significant
if maxDelta_x < SETTINGS.THRESHOLD
    
    % delta is insignificant. 
    t1 = upperLimit_k1;
    
    %if 
        % All subresultants are rank deficient (Singular)
    %    t1 = 0
    %else
        % All subresultants are full rank (non-singular)
    %    t1 = upperLimit_k1
    %end
    
    
else 
    t1 = lowerLimit_k1 + idx - 1;
end



delta_y = diff(log10(mat_metric),1,2);
vDelta_y = sum(delta_y,1);
[maxDelta_y, idx] = max(vDelta_y);

if maxDelta_y < SETTINGS.THRESHOLD
    
    % Delta is insignificant
    %t2 = 0;
    
    t2 = upperLimit_k2;
    
else
    % Delta is significant
    t2 = lowerLimit_k2 + idx - 1;
end



fprintf([mfilename ' : ' 'The Calculated Degree of the GCD is given by \n'])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt x : t1 = %i\n', t1)])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt y : t2 = %i\n', t2)])
end

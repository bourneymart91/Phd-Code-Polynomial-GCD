function [idx_optColumn] = GetOptimalColumn_2Polys(fxy, gxy, m, n, k, k1, k2)
% Get the optimal column c_{k} of the Sylvester subresultant matrix S(f,g)
% which when removed from S(f,g) gives minimal residual in the equation
% A_{k}x = c_{k} when solving for x.
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
% k : (Int) Total degree of polynomial d(x,y)
%
% k1 : (Int) Degree of polynomial d(x,y) with respect to x
%
% k2 : (Int) Degree of polynomial d(x,y) with respect to y
%
% % Outputs
%
% idx_col : (Int) Index of column to be removed from the Sylvester subresultant
% matrix.

global SETTINGS

switch SETTINGS.DEGREE_METHOD
    case 'Total'
        
        % Build the kth Sylvester matrix
        Sk = BuildT_Total_Bivariate_2Polys(fxy, gxy, m, n, k);
        
        % Get optimal column for removal
        idx_optColumn = GetOptimalColumn_Total(Sk);
        
    case 'Relative'
        
        % Build the kth Sylvester matrix
        Sk1k2 = BuildT_Relative_Bivariate_2Polys(fxy, gxy, k1, k2);
        
        % Get optimal column for removal
        idx_optColumn = GetOptimalColumn_Relative(Sk1k2);
        
    case 'Both'
        
        % Build the kth Sylvester matrix
        Skk1k2 = BuildT_Both_Bivariate_2Polys(fxy, gxy, m, n, k, k1, k2);
        
        % Get Optimal column for removal
        idx_optColumn = GetOptimalColumn_Both(Skk1k2);
        
    otherwise
        error('err')
end


end
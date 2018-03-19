function [idx_optColumn] = GetOptimalColumn_3Polys(fxy, gxy, hxy, m, n, o, k, k1, k2)
% Get the optimal column c_{k} of the Sylvester subresultant matrix S(f,g)
% which when removed from S(f,g) gives minimal residual in the equation
% A_{k}x = c_{k} when solving for x.
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomials f(x,y) 
%
% gxy : (Matrix) Coefficients of the polynomials g(x,y)
%
% hxy : (Matrix) Coefficients of the polynomials h(x,y)
%
% m n o : (Int) (Int) (Int) Total degree of f(x,y) g(x,y) and h(x,y)
% 
% k : (Int) Total degree of polynomial d(x,y)
%
% k1 : (Int) Degree of polynomial d(x,y) with respect to x
%
% k2 : (Int) Degree of polynomial d(x,y) with respect to y
%
% % Outputs
%
% idx_optColumn : Index of column to be removed from the Sylvester subresultant
% matrix.

global SETTINGS

switch SETTINGS.DEGREE_METHOD
    case 'Total'
        
        % Build the kth Sylvester matrix
        Sk = BuildT_Total_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, k);
        
        % Get optimal column for removal
        idx_optColumn = GetOptimalColumn_Total(Sk);
        
    case 'Relative'
        
        % Build the kth Sylvester matrix
        Sk1k2 = BuildT_Relative_Bivariate_3Polys(fxy, gxy, hxy, k1, k2);
        
        % Get optimal column for removal
        idx_optColumn = GetOptimalColumn_Relative(Sk1k2);
        
    case 'Both'
        
        % Build the kth Sylvester matrix
        Skk1k2 = BuildT_Both_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, k, k1, k2);
        
        % Get Optimal column for removal
        idx_optColumn = GetOptimalColumn_Both(Skk1k2);
        
    otherwise
        error('err')
end


end
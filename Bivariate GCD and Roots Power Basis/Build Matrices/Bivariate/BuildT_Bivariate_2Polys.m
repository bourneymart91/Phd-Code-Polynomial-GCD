function Sk = BuildT_Bivariate_2Polys(fxy, gxy, m, n, k, k1, k2)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomials f(x,y)  
%
% gxy : (Matrix) Coefficients of the polynomials g(x,y) 
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% k : (Int) Total degree of common divisor d(x,y)
%
% k1 : (Int) Degree of d(x,y) with respect to x
%
% k2 : (Int) Degree of d(x,y) with respect to y
%
% % Outputs
%
% Sk : (Matrix) kth Sylvester subresultant matrix



global SETTINGS


switch SETTINGS.DEGREE_METHOD
    
    case 'Total' 
        
        Sk = BuildT_Total_Bivariate_2Polys(fxy, gxy, m, n, k);
        
    case 'Relative'
        
        Sk = BuildT_Relative_Bivariate_2Polys(fxy, gxy, k1, k2);
        
    case 'Both'
        
        Sk = BuildT_Both_Bivariate_2Polys(fxy, gxy, m, n, k, k1, k2);
        
end

end
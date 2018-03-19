function Sk = BuildT_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, k, k1, k2)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
% 
% gxy : (Matrix) Coefficients of polynomial g(x,y) 
%
% hxy : (Matrix) Coefficients of polynomial h(x,y)
%
% m : (Int) Total degree of f(x,y) 
%
% n : (Int) Total degree of g(x,y)
%
% o : (Int) Total degree of h(x,y)
%
% k : (Int) Total degree of d(x,y)
%
% k1 : Relative degree of d(x,y) with respect to x
%
% k2 : Relative degree of d(x,y) with respect to y
%
% % Outputs
%
% Sk : (Matrix) kth Sylvester subresultant matrix of polynomials f, g and h

global SETTINGS


switch SETTINGS.DEGREE_METHOD
    
    case 'Total' 
        
        Sk = BuildT_Total_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, k);
        
    case 'Relative'
        
        Sk = BuildT_Relative_Bivariate_3Polys(fxy, gxy, hxy, k1, k2);
        
    case 'Both'
        
        Sk = BuildT_Both_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, k, k1, k2);
        
end

end
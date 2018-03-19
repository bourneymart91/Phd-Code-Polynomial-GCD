function [dxy] = GetGCDCoefficients_Bivariate_3Polys(fxy, gxy, hxy, uxy, vxy, wxy, m, n, o, k)
% GetGCDCoefficients(fxy,gxy,uxy,vxy,m,n,k)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of the polynomial h(x,y)
%
% uxy : (Matrix) Coefficients of the polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of the polynomial v(x,y)
%
% m n : (Int) (Int) : Total degree of polynomial f(x,y) and g(x,y)
%
% t : (Int) Total degree of polynomial d(x,y)

global SETTINGS

switch SETTINGS.DEGREE_METHOD
    case 'Total'
    
        dxy = GetGCDCoefficients_Total_Bivariate_3Polys(fxy, gxy, hxy, uxy, vxy, wxy, m, n, o, k);
        
    case 'Relative'
        
        dxy = GetGCDCoefficients_Relative_Bivariate_3Polys(fxy, gxy, hxy, uxy, vxy, wxy);
    
    case 'Both'
        
        dxy = GetGCDCoefficients_Both_Bivariate_3Polys(fxy, gxy, hxy, uxy, vxy, wxy, m, n, o, k);
        
    otherwise
        
        error([mfilename ' : DEGREE_METHOD must be Total, Relative, or Both'])
        
end
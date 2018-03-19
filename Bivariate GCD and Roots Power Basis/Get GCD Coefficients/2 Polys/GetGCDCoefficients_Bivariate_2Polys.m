function [dxy] = GetGCDCoefficients_Bivariate_2Polys(fxy, gxy, uxy, vxy, m, n, k)
% GetGCDCoefficients(fxy,gxy,uxy,vxy,m,n,k)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
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
    
        dxy = GetGCDCoefficients_Total_Bivariate_2Polys(fxy, gxy, uxy, vxy, m, n, k);
        
    case 'Relative'
        
        dxy = GetGCDCoefficients_Relative_Bivariate_2Polys(fxy, gxy, uxy, vxy);
    
    case 'Both'
        
        dxy = GetGCDCoefficients_Both_Bivariate_2Polys(fxy, gxy, uxy, vxy, m, n, k);
        
    otherwise
        
        error([mfilename ' : DEGREE_METHOD must be Total, Relative, or Both'])
        
end

end
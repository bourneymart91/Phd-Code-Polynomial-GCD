function [uxy, vxy, wxy] =  GetQuotients_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, t, t1, t2)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomials f(x,y)
%
% gxy : (Matrix) Coefficients of polynomials g(x,y)
%
% hxy : (Matrix) Coefficients of polynomials h(x,y)
%
% m n o : (Int) (Int) (Int) : Total degree of f(x,y), g(x,y) and h(x,y)
%
% t : Total degree of GCD d(x,y)
%
% t1 t2 : (Int) (Int) Relative degrees of d(x,y) with respect to x and y
%
% % Outputs
%
% uxy : (Matrix) Coefficients of polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of polynomial u(x,y)
%
% wxy : (Matrix) Coefficients of polynomial w(x,y)

global SETTINGS

% % Get coefficients of the quotients u(x,y) and v(x,y)
switch SETTINGS.DEGREE_METHOD
    
    case 'Total'
        
        [uxy, vxy, wxy] = GetQuotients_Total_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, t);
        
    case 'Relative'
        
        [uxy, vxy, wxy] = GetQuotients_Relative_Bivariate_3Polys(fxy, gxy, hxy, t1, t2);
        
    case 'Both'
        
        [uxy, vxy, wxy] = GetQuotients_Both_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, t, t1, t2);
        
    otherwise
        
        error('calc method is either total or relative')
        
end
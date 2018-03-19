function [uxy,vxy] =  GetQuotients_Bivariate_2Polys(fxy, gxy, m, n, t, t1, t2)
%
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomials g(x,y)
%
% m : (Int) Total degree of polynomial f(x,y)
%
% n : (Int) Total degree of polynomial g(x,y)
%
% t : (Int) Total degree of common divisor
%
% t1 t2 : (Int) (Int) : Relative degree of common divisor d(x,y)
% 
% % Outputs
%
% uxy : (Matrix) Coefficients of polynomials u(x,y)
%
% vxy : (Matrix) Coefficients of polynomials v(x,y)

global SETTINGS

% % Get coefficients of the quotients u(x,y) and v(x,y)
switch SETTINGS.DEGREE_METHOD
    case 'Total'
        
        [uxy,vxy] = GetQuotients_Total_Bivariate_2Polys(fxy, gxy, m, n, t);
        
    case 'Relative'
        
        [uxy,vxy] = GetQuotients_Relative_Bivariate_2Polys(fxy, gxy, t1, t2);
        
    case 'Both'
        
        [uxy,vxy] = GetQuotients_Both_Bivariate_2Polys(fxy, gxy, m, n, t, t1, t2);
        
    otherwise
        
        error([mfilename ': Calculation method must be total, relative or both'])
end



end
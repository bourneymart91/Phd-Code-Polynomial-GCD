function hxy = Deconvolve_Bivariate_Single(fxy, gxy, m, n)
% Performs a bivariate deconvolution of polynomials f(x,y) and g(x,y) to
% obtain h(x,y)
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

global SETTINGS
switch SETTINGS.DEGREE_METHOD
    case 'Total'
        
        hxy = Deconvolve_Bivariate_Single_Total(fxy,gxy,m,n);
        
    case 'Respective'
        
        hxy = Deconvolve_Bivariate_Single_Respective(fxy,gxy);
        
    case 'Both'
        
        hxy = Deconvolve_Bivariate_Single_Both(fxy,gxy,m,n);
        
    otherwise
        
        error('err')
end

end
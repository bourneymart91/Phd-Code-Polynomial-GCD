function [fxy_lr, gxy_lr, uxy_lr, vxy_lr] = STLN(fxy,gxy,m,n,k,k1,k2,idx_col)
% STLN(fxy,gxy,m,n,k,k1,k2,idx_col)
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) 
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% k : (Int) Total degree of d(x,y)
%
% k1 : (Int) Degree of d(x,y) with respect to x
%
% k2 : (Int) Degree of d(x,y) with respect to y
%
% % Outputs.
%
% fxy_lr : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy_lr : (Matrix) Coefficients of polynomial g(x,y)
%
% uxy_lr : (Matrix) Coefficients of polynomial u(x,y)
%
% vxy_lr : (Matrix) Coefficients of polynomial v(x,y)

global SETTINGS

switch SETTINGS.DEGREE_METHOD
    case 'Total'
        
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr] = STLN_Total(fxy, gxy, m, n, k, idx_col);
        
        
    case 'Relative'
        
        
        % Perform STLN to obtain low rank approximation
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr] = STLN_Relative(fxy, gxy, k1, k2, idx_col);
        
        
        %BuildT_Relative_Bivariate_2Polys(fxy, gxy, k1, k2);
        %BuildT_Relative_Bivariate_2Polys(fxy_lr, gxy_lr, k1, k2);
        
    case 'Both'
        
        % Perform STLN to obtain low rank approximation
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr] = STLN_Both(fxy, gxy, m, n, k, k1, k2, idx_col);
        
 
        
end


end
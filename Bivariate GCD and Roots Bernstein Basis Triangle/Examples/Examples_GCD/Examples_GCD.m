function [fxy,gxy,uxy,vxy,dxy,m,n,t] = Examples_GCD(ex_num)
%
% % Inputs
%
% ex_num : Example Number as a string
%
% % Outputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) 
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% uxy : (Matrix) Coefficients of polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of polynomial v(x,y)
%
% dxy : (Matrix) Coefficients of polynomial d(x,y)
%
% m : (Int) Degree of polynomial f(x,y)
%
% n : (Int) Degree of polynomial g(x,y)
%
% t : (Int) Degree of polynomial d(x,y)



[fxy, gxy, uxy, vxy, dxy, m, n, t] = ...
    Examples_GCD_FromCoefficients_2Polys(ex_num);
        
  

end

function [alpha,th1,th2] = Preprocess(fxy, gxy, m, n, k)
% Obtain optimal values of alpha and theta
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

% Get global variables
global SETTINGS

if (SETTINGS.BOOL_ALPHA_THETA)
    
    % Get maximum of each coefficient a_{i,j} from its entries in the Sylvester
    % Matrix S_{k}(f,g)
    [max_fxy, min_fxy] = GetMaxMin(fxy, m, n - k);
    
    % Get maximum of each coefficient b_{i,j} from its entries in the Sylvester
    % Matrix S_{k}(f,g)
    [max_gxy, min_gxy] = GetMaxMin(gxy, n, m-k);
    
    % Get Optimal alpha and theta
    [alpha, th1,th2] = OptimalAlphaTheta(max_fxy, min_fxy, max_gxy, min_gxy, m, n);
    
else
    alpha = 1;
    th1 = 1;
    th2 = 1;
    
end

end
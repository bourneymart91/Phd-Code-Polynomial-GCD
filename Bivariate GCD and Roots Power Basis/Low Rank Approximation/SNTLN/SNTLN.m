function [fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr,th1_lr,th2_lr] = ...
    SNTLN(fxy, gxy, alpha, th1, th2, m, n, t, t1, t2, idx_col)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% alpha : (Float) Optimal value of \alpha
%
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}
%
% m : (Int) Total degree of polynomial f(x,y)
%
% n : (Int) Total degree of polynomial g(x,y)
%
% t : (Int) Total degree of polynomial d(x,y)
%
% t1 : (Int) Degree of d(x,y) with respect to x
%
% t2 : (Int) Degree of d(x,y) with respcet to y
%
% idx_col : (Int) Index of the optimal column to be removed from the Sylvester
% subresultant matrix.
%
% % Outputs.
%
% fxy_lr : (Matrix) Coefficients of polynomial f(x,y) + \delta f(x,y)
%
% gxy_lr : (Matrix) Coefficients of polynomial g(x,y) + \delta g(x,y)
%
% uxy_lr : (Matrix) Coefficients of polynomial u(x,y) + \delta u(x,y)
%
% vxy_lr : (Matrix) Coefficients of polynomial g(x,y) + \delta g(x,y)
%
% alpha_lr : (Float) \alpha + \delta \alpha
%
% th1_lr : (Float) \theta_{1} + \delta \theta_{1}
%
% th1_lr : (Float) \theta_{2} + \delta \theta_{2}
%
% x_lr : (Vector) The solution vector x in the problem Ax = b.

% Initialise Global Settings
global SETTINGS

% Get method to compute SNTLN
switch SETTINGS.DEGREE_METHOD
    case 'Total'
        %
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
            SNTLN_Total(fxy, gxy, m, n, alpha, th1, th2, t, idx_col);
        
    case 'Relative'
        % Get the SNTLN of the Sylvester matrix
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
            SNTLN_Relative(fxy, gxy, alpha, th1, th2, t1, t2, idx_col);
        
    case 'Both'
        
        % Get the SNTLN of the Sylvester matrix
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
            SNTLN_Both(fxy, gxy, alpha, th1, th2, m, n, t, t1, t2, idx_col);
end



end







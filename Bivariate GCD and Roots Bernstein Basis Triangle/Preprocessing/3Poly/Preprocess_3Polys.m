function [lambda, mu, rho, theta1, theta2] = Preprocess_3Polys(fxy, gxy, ...
    hxy, m, n, o, k)




global SETTINGS

switch SETTINGS.SYLVESTER_MATRIX_3POLY_N_EQUATIONS
    
    
    case '2'
        [lambda, mu, rho, theta1, theta2] = Preprocess_3Polys_2Eqns(fxy, gxy, ...
            hxy, m, n, o, k);
        
    case '3'
        
        [lambda, mu, rho, theta1, theta2] = Preprocess_3Polys_3Eqns(fxy, gxy, ...
            hxy, m, n, o, k);
        
end

end





function [lambda, mu, rho, theta1, theta2] = Preprocess_3Polys_2Eqns(fxy, gxy, ...
    hxy, m, n, o, k)
% Obtain optimal values of alpha and theta
%
% % Inputs.
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
% % Outputs
%
% alpha : (Float)
%
% beta : (Float)
%
% gamma : (Float)
%
% theta1 : (Float)
%
% theta2 : (Float)

% Get global variables
global SETTINGS




if (SETTINGS.BOOL_ALPHA_THETA)
    
    % Get maximum of each coefficient a_{i,j} from its entries in the Sylvester
    % Matrix S_{k}(f,g)
    [max_fxy1, min_fxy1] = GetMaxMin(fxy, m, n - k);
    [max_fxy2, min_fxy2] = GetMaxMin(fxy, m, o - k);
    [max_gxy, min_gxy] = GetMaxMin(gxy, n, m - k);
    [max_hxy, min_hxy] = GetMaxMin(hxy, o, m - k);
    
    
    
    
    % Get Optimal alpha and theta
    [lambda, mu, rho, theta1, theta2] = ...
        OptimalAlphaBetaTheta(max_fxy1, min_fxy1,...
        max_fxy2, min_fxy2, ...
        max_gxy, min_gxy, ...
        max_hxy, min_hxy,...
        m, n, o);
    
else
    lambda = 1;
    mu = 1;
    rho = 1;
    theta1 = 1;
    theta2 = 1;
    
end

end



function [lambda, mu, rho, theta1, theta2] = Preprocess_3Polys_3Eqns(fxy, gxy, ...
    hxy, m, n, o, k)
% Obtain optimal values of alpha and theta
%
% % Inputs.
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
% % Outputs
%
% alpha : (Float)
%
% beta : (Float)
%
% gamma : (Float)
%
% theta1 : (Float)
%
% theta2 : (Float)

% Get global variables
global SETTINGS




if (SETTINGS.BOOL_ALPHA_THETA)
    
    % Get maximum of each coefficient a_{i,j} from its entries in the Sylvester
    % Matrix S_{k}(f,g)
    [max_fxy1, min_fxy1] = GetMaxMin(fxy, m, n - k);
    [max_fxy2, min_fxy2] = GetMaxMin(fxy, m, o - k);
    [max_gxy1, min_gxy1] = GetMaxMin(gxy, n, m - k);
    [max_gxy2, min_gxy2] = GetMaxMin(gxy, m, o - k);
    [max_hxy1, min_hxy1] = GetMaxMin(hxy, o, m - k);
    [max_hxy2, min_hxy2] = GetMaxMin(hxy, o, n - k);
    
    
    
    % Get Optimal alpha and theta
    [lambda, mu, rho, theta1, theta2] = ...
        OptimalAlphaBetaTheta_3Polys_3Eqns(...
        max_fxy1, min_fxy1, ...
        max_fxy2, min_fxy2, ...
        max_gxy1, min_gxy1, ...
        max_gxy2, min_gxy2, ...
        max_hxy1, min_hxy1, ...
        max_hxy2, min_hxy2, ...
        m, n, o);
    
else
    lambda = 1;
    mu = 1;
    rho = 1;
    theta1 = 1;
    theta2 = 1;
    
end

end
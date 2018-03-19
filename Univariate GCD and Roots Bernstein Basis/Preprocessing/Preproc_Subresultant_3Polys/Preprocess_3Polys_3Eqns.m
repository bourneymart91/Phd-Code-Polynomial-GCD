function [GM_fx, GM_gx, GM_hx, alpha, beta, gamma, theta] = ...
    Preprocess_3Polys_3Eqns(fx, gx, hx, k)
% Preprocess_3Polys(fx, gx, hx, k)
%
% Get the optimal values of lamdba, mu, alpha and theta.
%
% % Inputs.
%
%
% fx : (Vector) Vector of coefficients of f(x)
%
% gx : (Vector) Vector of coefficients of g(x)
%
% hx : (Vector) Vector of coefficients of h(x)
%
% k : (Int) Degree of common divisor d_{k}(x)
%
% % Outputs.
%
% GM_fx : (Float) Geometric mean of entries of f(x) in k-th subresultant matrix
%
% GM_gx : (Float) Geometric mean of entries of g(x) in k-th subresultant matrix
%
% GM_hx : (Float) Geometric mean of entries of h(x) in k-th subresultant matrix
%
% alpha : (Float) Optimal value of \alpha
%
% beta : (Float) Optimal value of \beta
%
% theta : (Float) Optimal value of \theta


% Global variables
global SETTINGS

% Get degree of polynomials f(x), g(x) and h(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

% Get the mean of the entries of f(x) in T_{n-k}(f) and T_{o-k}(f)


GM_fx = GetMean_3Polys_3Eqns(fx, n - k, o - k);
GM_gx = GetMean_3Polys_3Eqns(gx, m - k, o - k);
GM_hx = GetMean_3Polys_3Eqns(hx, m - k, n - k);

%GM_fx_test1 = GetMean(fx, n - k);
%GM_fx_test2 = GetMean_3PolySubresultant(fx, n - k, o - k);

% Normalize f(x) and g(x) by geometric means
fx_n = fx ./ GM_fx;
gx_n = gx ./ GM_gx;
hx_n = hx ./ GM_hx;


if (SETTINGS.BOOL_ALPHA_THETA)
    
    % For each coefficient ai of F, obtain the max and min such that F_max =
    % [max a0, max a1,...] and similarly for F_min, G_max, G_min
    
    % Get maximum and minimum entries of each a_{i} in the first
    % partition T_{n-k}(f)
    [v_F_max1,v_F_min1] = GetMaxMin(fx_n, n - k);
    [v_F_max2, v_F_min2] = GetMaxMin(fx_n, o - k);
    [v_G_max1, v_G_min1] = GetMaxMin(gx_n, m - k);
    [v_G_max2, v_G_min2] = GetMaxMin(gx_n, o - k);
    [v_H_max1, v_H_min1] = GetMaxMin(hx_n, m - k);
    [v_H_max2, v_H_min2] = GetMaxMin(hx_n, n - k);
    
    % Calculate the optimal value of alpha and theta for the kth
    % subresultant matrix.

    [alpha, beta, gamma, theta] = ...
        OptimalAlphaBetaGammaTheta_3Polys_3Eqns(...
        v_F_max1, v_F_min1, ...
        v_F_max2, v_F_min2, ...
        v_G_max1, v_G_min1, ...
        v_G_max2, v_G_min2, ...
        v_H_max1, v_H_min1, ...
        v_H_max2, v_H_min2);
   
  
    
else
    alpha = 1;
    beta = 1;
    gamma = 1;
    theta = 1;
    
end

end

function [] = PrintToFile(F_max,F_min,G_max,G_min,m,n,k,alpha,theta)

global SETTINGS

fullFileName = 'Results_Preprocessing.txt';


if exist('Results_Preprocessing.txt', 'file')
    fileID = fopen('Results_Preprocessing.txt','a');
    fprintf(fileID,'%s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s,\t %s,\t %s, \t %s \n',...
        datetime('now'),...
        SETTINGS.PROBLEM_TYPE,...
        SETTINGS.EX_NUM,...
        int2str(m),...
        int2str(n),...
        int2str(k),...
        SETTINGS.MEAN_METHOD,...
        F_max,...
        F_min,...
        G_max,...
        G_min,...
        num2str(alpha),...
        num2str(theta)...
        );
    fclose(fileID);
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end

end

function [] = print(v_F_max,v_F_min,v_G_max,v_G_min,m,n,k)
% Get maximum entry of all entries of T_{n-k}(f)
f_max = max(v_F_max);

% Get minimum entry of all entries of T_{n-k}(f)
f_min = min(v_F_min);

% Get maximum entry of all entries of T_{m-k}(g)
g_max = max(v_G_max);

% Get minimum entry of all entries of T_{m-k}(g)
g_min = min(v_G_min);

% Print max and minimum entries
PrintToFile(f_max,f_min,g_max,g_min,m,n,k,1,1)
end
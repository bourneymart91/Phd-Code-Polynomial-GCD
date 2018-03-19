function [GM_fx, GM_gx, alpha, theta] = Preprocess_2Polys(fx, gx, k)
% Preprocess_2Polys(fx, gx, k)
% Get the optimal values of lamdba, mu, alpha and theta.
%
% Inputs.
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% gx : (Vector) Coefficients of polynomial g(x)
%
% k : (Int) Degree of common divisor d_{k}(x)
%
% Outputs.
%
% GM_fx : (Float) Geometric mean of entries of f(x) in k-th subresultant matrix
%
% GM_gx : (Float) Geometric mean of entries of g(x) in k-th subresultant matrix
%
% alpha : (Float) Optimal value of alpha
%
% theta : (Float) Optimal value of theta


% Global variables
global SETTINGS

% Get degree of polynomial f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Get the mean of the entries of T_{n-k}(F) and T_{m-k}(g)
GM_fx = GetMean(fx, n - k);
GM_gx = GetMean(gx, m - k);

% Normalize f(x) and g(x) by geometric means
fx_n = fx./ GM_fx;
gx_n = gx./ GM_gx;





if (SETTINGS.BOOL_ALPHA_THETA)
    
    
    % For each coefficient ai of F, obtain the max and min such that F_max =
    % [max a0, max a1,...] and similarly for F_min, G_max, G_min
    
    
    [v_F_max, v_F_min] = GetMaxMin(fx_n, n - k);
    [v_G_max, v_G_min] = GetMaxMin(gx_n, m - k);
    
    
    % Calculate the optimal value of alpha and theta for the kth
    % subresultant matrix.
    [alpha, theta] = OptimalAlphaTheta(v_F_max, v_F_min, v_G_max, v_G_min);
  
   
else
    alpha = 1;
    theta = 1;
end

end

function [] = PrintToFile(F_max, F_min, G_max, G_min, m, n, k, alpha, theta)

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

function [] = print(v_F_max, v_F_min, v_G_max, v_G_min, m, n, k)

% Get maximum entry of all entries of T_{n-k}(f)
f_max = max(v_F_max);

% Get minimum entry of all entries of T_{n-k}(f)
f_min = min(v_F_min);

% Get maximum entry of all entries of T_{m-k}(g)
g_max = max(v_G_max);

% Get minimum entry of all entries of T_{m-k}(g)
g_min = min(v_G_min);

% Print max and minimum entries
PrintToFile(f_max, f_min, g_max, g_min, m, n, k, 1, 1)
end
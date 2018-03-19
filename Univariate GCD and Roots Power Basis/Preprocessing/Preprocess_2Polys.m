function [GM_fx, GM_gx, alpha, theta] = Preprocess_2Polys(fx, gx)
% Preprocess(fx,gx)
%
% Preprocess the input polynomials to obtain geometric means of the
% coefficients of f(x) and g(x), and optimal values of alpha and theta such
% that the ratio of entries in sylvester matrix S(f,g) is minimized.
%
% Inputs.
%
% [fx, gx] : Vector of coefficients of polynomial f(x) and g(x)
%
% % Outputs
%
% [GM_fx, GM_gx]
%
% alpha
%
% theta

global SETTINGS

m = GetDegree(fx);
n = GetDegree(gx);


% Get Mean of entries of f(x) in C_{n-k}(f)
GM_fx = GetMean(fx);

% Get Mean of entries of g(x) in C_{m-k}(g)
GM_gx = GetMean(gx);

% Divide f(x) and g(x) by geometric mean
fx_n = fx./ GM_fx;
gx_n = gx./ GM_gx;

if(SETTINGS.BOOL_ALPHA_THETA)
    
    
    % Get opitmal values of alpha and theta
    [alpha, theta] = GetOptimalAlphaAndTheta_2Polys(fx_n, gx_n);
    
    % Obtain f(w) and g(w) from f(x) and g(x)]
    fw = GetWithThetas(fx_n, theta);
    gw = GetWithThetas(gx_n, theta);
    
    F_max = max(fx_n);
    F_min = min(fx_n);
    G_max = max(gx_n);
    G_min = min(gx_n);
    PrintToFile(F_max, F_min, G_max, G_min, m, n, '1', '1');
    
    %%
    F_max = max(fw);
    F_min = min(fw);
    G_max = max(gw);
    G_min = min(gw);
    
    PrintToFile(F_max,F_min,G_max,G_min,m,n,alpha,theta);
    

else
    
    % Dont apply preprocessing
    theta = 1;
    alpha = 1;
    
    
end

end



function [] = PrintToFile(F_max,F_min,G_max,G_min,m,n,alpha,theta)

global SETTINGS

fullFileName = 'Preprocessing/Results_Preprocessing.txt';


if exist(fullFileName, 'file')
    fileID = fopen(fullFileName,'a');
    fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
        datetime('now'),...
        SETTINGS.EX_NUM,...
        int2str(m),...
        int2str(n),...
        SETTINGS.MEAN_METHOD,...
        F_max,...
        F_min,...
        G_max,...
        G_min,...
        SETTINGS.BOOL_ALPHA_THETA,...
        alpha,...
        theta...
        );
    fclose(fileID);
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end

end
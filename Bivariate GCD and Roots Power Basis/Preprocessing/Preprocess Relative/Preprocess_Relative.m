function [lambda,mu,alpha,th1,th2] = Preprocess_Relative(fxy, gxy)
% Preprocess the polynomial f(x,y) and g(x,y), and return the geometric
% means and optimal alphas and thetas from the preprocessing.
%
% Inputs.
%
% fxy : Matrix of coefficients of polynomial f(x,y)
%
% gxy : Matrix of coefficients of polynomial g(x,y)
%
%
% Outputs.
%
%
% lambda : Geometric mean of coefficients of f(x,y)
%
% mu : Geometric mean of coefficients of g(x,y)
%
% alpha : Optimal \alpha
%
% th1 : Optimal \theta_{1}
%
% th2 : Optimal \theta_{2}

% Global variables

global SETTINGS

% Get degree of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degree of g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Get mean of f(x,y) in entries of C_{}(f)
lambda = GetMean(fxy);

% Get mean of g(x,y) in entries of C_{}(g)
mu = GetMean(gxy);

% Normalise coefficients of f(x,y) and g(x,y) by dividing the
% coefficients by the respective means.
fxy_n = fxy ./ lambda;
gxy_n = gxy ./ mu;

switch SETTINGS.BOOL_ALPHA_THETA
    case true
        
        
        optimization_method = 'Together';
        
        %optimization_method = 'Independent';
        
        switch optimization_method
            case 'Independent'
                
                % Obtain optimal values of alpha, theta_{1} and theta_{2}
                alpha = OptimalAlpha(fxy_n,gxy_n);
                [th1,th2] = OptimalTheta(fxy_n,alpha.*gxy_n);
                
            case 'Together'
                
                [alpha, th1,th2] = OptimalAlphaAndTheta_Relative(fxy_n,gxy_n);
                
                str1 = sprintf('Alpha : %2.4f', alpha);
                str2 = sprintf('Theta_{1} : %2.4f', th1);
                str3 = sprintf('Theta_{2} : %2.4f', th2);
                
                fprintf([mfilename ' : ' str1 '\n']);
                fprintf([mfilename ' : ' str2 '\n']);
                fprintf([mfilename ' : ' str3 '\n']);
                
                
        end
        
        
        % Get f(w,w) from f(x,y)
        fww = GetWithThetas(fxy_n,th1,th2);
        
        % Get g(w,w) from g(x,y)
        gww = GetWithThetas(gxy_n,th1,th2);
        
        
        
        % Get the coefficients of f(x,y) and f(w,w) as vectors
        v_fxy = GetAsVector(fxy_n);
        v_fww = GetAsVector(fww);
        
        % Get the coefficients of g(x,y) and g(w,w) as vectors
        v_gxy = GetAsVector(gxy_n);
        v_gww = GetAsVector(gww);
        
        
        
        if (SETTINGS.PLOT_GRAPHS)
            % Plot the coefficients of f(x,y) and f(w,w)
            PlotCoefficients(v_fxy,v_fww,'f')
            % Plot the coefficients of g(x,y) and g(w,w)
            PlotCoefficients(v_gxy,alpha.*v_gww,'g')
        end
        
        % Get maximum coefficient of f(x)
        max_fx = max(abs(v_fxy(v_fxy~=0)));
        % Get minimum coefficient of f(x)
        min_fx = min(abs(v_fxy(v_fxy~=0)));
        
        % Get maximum coefficient of g(x)
        max_gx = max(abs(v_gxy(v_gxy~=0)));
        % Get minimum coefficient of g(x)
        min_gx = min(abs(v_gxy(v_gxy~=0)));
        
        PrintToFile([m1,m2], [n1,n2], max_fx, min_fx, max_gx, min_gx, 1, 1, 1);
        
        
        % Get maximum coefficient of f(x)
        max_fw = max(abs(v_fww(v_fww~=0)));
        % Get minimum coefficient of f(x)
        min_fw = min(abs(v_fww(v_fww~=0)));
        
        
        
        a_gw = alpha.*v_gww;
        
        % Get maximum coefficient of g(x)
        max_gw = max(abs(a_gw(a_gw~=0)));
        
        % Get minimum coefficient of g(x)
        min_gw = min(abs(a_gw(a_gw~=0)));
        
        PrintToFile([m1,m2],[n1,n2],max_fw,min_fw,max_gw,min_gw,alpha,th1,th2);
        
        
        str1 = sprintf('Condition S(f(x),g(x)) : %2.4e', cond(BuildT_Relative_Bivariate_2Polys(fxy, gxy, 1, 1)));
        str2 = sprintf('Condition S(f(w),alpha*g(w)) : %2.4e', cond(BuildT_Relative_Bivariate_2Polys(fww, alpha.*gww, 1, 1)));
        
        fprintf([mfilename ' : ' str1 '\n']);
        fprintf([mfilename ' : ' str2 '\n']);
        
        
    case false
        
        % Set linprog outputs to be 1
        th1 = 1;
        th2 = 1;
        alpha = 1;
        
    otherwise
        error('BOOL_ALPHA_THETA is either true or false')
end
end




function PrintToFile(pair_m, pair_n, max_fx, min_fx, max_gx, min_gx,alpha,th1,th2)
global SETTINGS

m1 = pair_m(1);
m2 = pair_m(2);
n1 = pair_n(1);
n2 = pair_n(2);

fullFileName = 'Results/Results_Preproc.txt';


if exist(fullFileName, 'file')
    fileID = fopen(fullFileName,'a');
    WriteNewLine()
    fclose(fileID);
else
    
    % File does not exist.
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    WriteNewLine()
    fclose(fileID);
    
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end

    function WriteNewLine()
        % Write a new line to the results file
        
        % 13 fields
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s  \n',...
            datetime('now'),...
            SETTINGS.EX_NUM,...
            int2str(m1),...
            int2str(m2),...
            int2str(n1),...
            int2str(n2),...
            max_fx,...
            min_fx,...
            max_gx,...
            min_gx,...
            num2str(alpha),...
            num2str(th1),...
            num2str(th2)...
            );
    end

    function WriteHeader()
        fprintf(fileID,'DATE,EX_NUM,m1,m2,n1,n2,max_fx,min_fx,max_gx,min_gx,alpha,th1,th2');
    end
end

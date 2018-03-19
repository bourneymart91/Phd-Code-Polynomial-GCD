function [] = o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, rank_revealing_metric)
% o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, rank_revealing_metric)
%
% Given two polynomials f(x) and g(x) calculate the GCD d(x).
%
% % Inputs.
%
% ex_num : (String) Example Number
%
% el : (Float) Noise lower level
%
% em : (Float) Noise upper level
%
% mean_method : (String)
%       'None'
%       'Geometric Mean Matlab Method'
%       'Geometric Mean My Method'
%
% bool_alpha_theta (Boolean)
%       * true :
%       * false :
%
% low_rank_approx_method: (String) 
%       'Standard SNTLN'
%       'Standard STLN'
%       'None'
%
% apf_method :(String)
%       'Standard Nonlinear'
%       'Standard Linear'
%       'None'
%
% rank_revealing_metric
%   * 'Minimum Singular Values'
%   * 'Residuals'
%   * 'R1 Row Diagonals'
%   * 'R1 Row Norms'
%
%
%
%
% % Example
%
% >> o_gcd_Univariate_2Polys('1', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'None', 'None', 'Minimum Singular Values')
% >> o_gcd_Univariate_2Polys('1', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'Standard STLN', 'Standard APF Nonlinear', 'Minimum Singular Values')

% Initialise global variables
global SETTINGS



% Add subfolders
restoredefaultpath

% Add that folder plus all subfolders to the path.
addpath(genpath(pwd));

% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end

% Set global variables
SetGlobalVariables_GCDFinding(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, rank_revealing_metric)

% Get coefficients of f(x,y) g(x,y) from example file
[fx_exact, gx_exact, dx_exact,ux_exact, vx_exact] = Examples_GCD(ex_num);

% Add noise to the coefficients of polynomials f(x) and g(x) at a
% predefined signal-to-noise ratio.
[fx_noisy, ~] = AddVariableNoiseToPoly(fx_exact, emin, emax);
[gx_noisy, ~] = AddVariableNoiseToPoly(gx_exact, emin, emax);

% Get degree of polynomials f(x,y) and g(x,y)
m = GetDegree(fx_noisy);
n = GetDegree(gx_noisy);

LineBreakLarge()
fprintf([mfilename ' : ' sprintf('Degree of the GCD is : %s \n', int2str(GetDegree(dx_exact)))])
LineBreakLarge()

fx = fx_noisy;
gx = gx_noisy;

% %
% %
% %
% Get the GCD by Zeng method
%[d_Zeng,v_Zeng,w_Zeng,res,cond] = o_gcd_zeng(fx,gx);

% %
% %
% %
% Get the GCD d(x) of f(x) and g(x) by my method

% Get upper and lower bound of degree of GCD.
lowerLimit_t = 1;
upperLimit_t = min(m, n);
limits_t = [lowerLimit_t, upperLimit_t];

rank_range = [0 0];

% Compute degree of gcd by my method
[fx_calc, gx_calc, dx_calc, ux_calc, vx_calc, ~, ~] = ...
    o_gcd_mymethod_Univariate_2Polys(fx, gx, limits_t, rank_range);

% Get distance of the computed d(x) from the exact d(x)

my_error.dx = GetDistance(dx_exact, dx_calc);
my_error.ux = GetDistance(ux_exact, ux_calc);
my_error.vx = GetDistance(vx_exact, vx_calc);

my_error.MyMethod = my_error.dx;

fprintf([mfilename sprintf(': Error d(x) : %2.4e \n', my_error.dx)])
fprintf([mfilename sprintf(': Error u(x) : %2.4e \n', my_error.ux)])
fprintf([mfilename sprintf(': Error v(x) : %2.4e \n', my_error.vx)])
    


PrintToFile(GetDegree(fx), GetDegree(gx), GetDegree(dx_exact), GetDegree(dx_calc), my_error)

% %
% %
% %
% Get the GCD by an alternative method
%[dx] = o_gcd_experiment_method(fx,gx)


% % Plot the three curves f(x), g(x) and d(x)

if( SETTINGS.PLOT_GRAPHS)
    
        % plot f(x) and g(x)
        t = linspace(-10,10,200);
        f_y = polyval(fx,t);
        g_y = polyval(gx,t);
        d_y = polyval(dx_calc,t);
        figure('name','Curve Plot')
        hold on
        plot(t,f_y,'DisplayName','f(y)')
        plot(t,g_y,'DisplayName','g(y)')
        plot(t,d_y,'DisplayName','d(y)')
        legend(gca,'show');
        hold off
    
end




end


function [dist] = GetDistance(f_exact, f_computed)
% GetDistance(u_exact,u_computed,name)
%
% Get the distance between the coefficients of two vectors.
%
% Inputs.
%
% f_exact : (Vector) Coefficients of f(x) exactly
%
% f_computed : (Vector) Coefficients of f(x) as computed
%
% name : (String) Name of function used for printing

% Normalise f(x) and f(x) computed.
f_exact = Normalise(f_exact);
f_computed = Normalise(f_computed);

% Get distance

try
    
    dist = norm(f_exact - f_computed) ./ norm(f_exact);
    
catch
    dist = 1000;
end
end


function [] = PrintToFile(m, n, t_exact, t_comp, error)
%
%
% % Inputs
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% t : Computed degree of the GCD d(x)
%
% error : array of errors e
%   error.dx
%   error.ux
%   error.vx

global SETTINGS

fullFileName = sprintf('Results/Results_o_gcd.txt');

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    WriteNewLine()
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    WriteNewLine()
    fclose(fileID);
    
end

    function WriteNewLine()
        
        % 19 fields
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
            datetime('now'),...,...
            SETTINGS.EX_NUM,...
            int2str(m),...
            int2str(n),...
            int2str(t_exact),...    
            int2str(t_comp),...
            num2str(error.ux),...
            num2str(error.vx),...
            num2str(error.dx),...
            SETTINGS.MEAN_METHOD,...
            num2str(SETTINGS.BOOL_ALPHA_THETA),...
            SETTINGS.EMIN,...
            SETTINGS.EMAX,...
            SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
            num2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
            SETTINGS.APF_METHOD,...
            num2str(SETTINGS.APF_REQ_ITE),...
            SETTINGS.GCD_COEFFICIENT_METHOD...
            );
            
    end

    function WriteHeader()
        fprintf(fileID,'DATE,EX_NUM,m,n,t_exact,t_comp,ERROR_UX,ERROR_VX,ERROR_DX,MEAN_METHOD,BOOL_ALPHA_THETA, EMIN, EMAX, LOW_RANK_APPROX_METHOD,LOW_RANK_ITE, APF_METHOD, APF_ITE,GCD_METHOD \n');
    end







end


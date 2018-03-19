function [] = o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, sylvester_build_method, low_rank_approx_method, ...
    apf_method, rank_revealing_metric, nEquations)
% o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, sylvester_build_method, low_rank_approx_method, apf_method, rank_revealing_metric)
%
% Given THREE polynomials f(x) and g(x) and h(x) calculate the GCD d(x).
%
% % Inputs.
%
% ex_num : (String) Example Number
%
% el : Noise lower level
%
% em : Noise upper level
%
% mean_method :
%       'None'
%       'Geometric Mean Matlab Method'
%       'Geometric Mean My Method'
%
% bool_alpha_theta : true false
%       true  : Include preproc
%   `   false : Exclude preproc
%
% sylvester_build_method
%       'T'
%       'T All Equations'
%
% low_rank_approx_method:
%       'Standard SNTLN'
%       'Standard STLN'
%       'None'
%
% apf_method :
%       'Standard Nonlinear'
%       'Standard Linear'
%       'None'
%
% rank_revealing_metric : 
%       'Minimum Singular Values'
%
% deconvolution_method : 
%
%
%
%
%
%
% % Example
% >> o_gcd_Univariate_3Polys('1', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'T All Equations', 'None', 'None', 'Minimum Singular Values'  )
% >> o_gcd_Univariate_3Polys('1', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'T All Equations', 'Standard STLN', 'Standard APF Nonlinear', 'Minimum Singular Values' )



% % Add subfolders
restoredefaultpath

addpath(genpath(pwd));

% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end

% Set global variables
SetGlobalVariables_GCD_3Polys(ex_num, emin, emax, ...
    mean_method, bool_alpha_theta, sylvester_build_method, low_rank_approx_method, apf_method, ...
    rank_revealing_metric, nEquations)

% Get coefficients of f(x,y) g(x,y) from example file
[fx_exact, gx_exact, hx_exact, dx_exact, ux_exact, vx_exact, wx_exact] = ...
    Examples_GCD_FromCoefficients_3Polys(ex_num);



% Add noise to the coefficients of polynomials f(x) and g(x) at a
% predefined signal-to-noise ratio.
[fx_noisy,~] = AddVariableNoiseToPoly(fx_exact, emin, emax);
[gx_noisy,~] = AddVariableNoiseToPoly(gx_exact, emin, emax);
[hx_noisy,~] = AddVariableNoiseToPoly(hx_exact, emin, emax);

%fx_noisy = fx_noisy ./ norm(fx_noisy);
%gx_noisy = gx_noisy ./ norm(gx_noisy);
%hx_noisy = hx_noisy ./norm(hx_noisy);




fx = fx_noisy ;
gx = gx_noisy ;
hx = hx_noisy ;

% %
% %
% %
% Get the GCD d(x) of f(x) and g(x) by my method

% Get upper and lower bound of degree of GCD.
upper_bound = min([GetDegree(fx), GetDegree(gx), GetDegree(hx)]);
lower_bound = 1;
limits_t = [lower_bound, upper_bound];

rank_range = [0 0];

% Compute degree of gcd by my method
[fx_calc, gx_calc, hx_calc, ...
    dx_calc, ...
    ux_calc, vx_calc, wx_calc, ~, ~] = o_gcd_mymethod_Univariate_3Polys(fx, gx, hx, limits_t, rank_range);

% Get distance of the computed d(x) from the exact d(x)

my_error.dx = GetDistance(dx_exact, dx_calc);
my_error.ux = GetDistance(ux_exact, ux_calc);
my_error.vx = GetDistance(vx_exact, vx_calc);
my_error.wx = GetDistance(wx_exact, wx_calc);

my_error.MyMethod = my_error.dx;

fprintf([mfilename sprintf(': Error d(x) : %2.4e \n', my_error.dx)])
fprintf([mfilename sprintf(': Error u(x) : %e \n', my_error.ux)])
fprintf([mfilename sprintf(': Error v(x) : %e \n', my_error.vx)])
fprintf([mfilename sprintf(': Error w(x) : %e \n', my_error.wx)])

PrintToFile(GetDegree(fx),GetDegree(gx),GetDegree(hx),GetDegree(dx_exact),GetDegree(dx_calc),my_error)

% %
% %
% %
% Get the GCD by an alternative method
%[dx] = o_gcd_experiment_method(fx,gx)




end


function [dist] = GetDistance(f_exact,f_computed)
% GetDistance(u_exact,u_computed,name)
%
% Get the distance between the coefficients of two vectors.
%
% Inputs.
%
% f_exact :
%
% f_computed :
%
% name : Name of function used for printing

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


function [] = PrintToFile(m,n,o,t_exact,t_comp,error)
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

fullFileName = sprintf('Results/Results_o_gcd_3Polys%s.txt',datetime('today'));

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
        fprintf(fileID,'%s,%s,%s, %s,%s,%s, %s,%s,%s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s \n',...
            datetime('now'),...,...
            SETTINGS.EX_NUM,...
            int2str(m),...
            int2str(n),...
            int2str(o),...
            int2str(t_exact),...    
            int2str(t_comp),...
            num2str(error.ux),...
            num2str(error.vx),...
            num2str(error.wx),...
            num2str(error.dx),...
            SETTINGS.MEAN_METHOD,...
            SETTINGS.BOOL_ALPHA_THETA,...
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


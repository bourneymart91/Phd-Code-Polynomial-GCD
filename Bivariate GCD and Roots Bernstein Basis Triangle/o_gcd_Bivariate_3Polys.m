function [] = o_gcd_Bivariate_3Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_type,...
    rank_revealing_metric, nEquations)
% o_gcd(ex_num, emin, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_type, rank_revealing_metric)
%
% Calculate the GCD d(x,y) of two polynomials f(x,y) and g(x,y) taken from
% the example file.
%
% % Inputs.
%
% ex_num  : (String) Example Number (String)
%
% emin : (Float) Minimum Noise level
%
% emax : (Float) Maximum signal to noise ratio
%
% mean_method : (String)
%       'Geometric Mean Matlab Method'
%       'None'
%
% bool_alpha_theta (Boolean)
%       true - Include Preprocessing
%       false - Exclude Preprocessing
%
% low_rank_approx_method (String)
%       'Standard SNTLN'
%       'None'
%
% apf_method (String)
%       'None'
%       'Standard APF Linear'
%       'Standard APF Nonlinear'
%
%
% sylvester_matrix_type : (String)
%       * T
%       * DT
%       * DTQ
%       * TQ
%
% rank_revealing_metric : (String)
%   'Minimum Singular Values'
%
%
% % Example
%
% >> o_gcd_Bivariate_3Polys('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'None', 'None', 'DTQ')
% >> o_gcd_Bivariate_3Polys('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'Standard STLN', 'Standard Nonlinear APF', 'DTQ')


% add path
restoredefaultpath
addpath(genpath(pwd));

% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end

% Set global variables
SetGlobalVariables_GCD_3Polys( ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_type, rank_revealing_metric, nEquations);

% Print Parameters to screen
fprintf('INPUTS. \n')
fprintf('EXAMPLE NUMBER %s \n',ex_num)
fprintf('EMIN : %s \n',emin)
fprintf('EMAX : %s \n',emax)
fprintf('MEAN METHOD : %s \n', mean_method)
fprintf('PREPROCESSING : %s \n',num2str(bool_alpha_theta))
fprintf('LOW RANK METHOD : %s \n',low_rank_approx_method)
fprintf('APF METHOD : %s \n', apf_method)
fprintf('SYLVESTER TYPE : %s \n', sylvester_type)
fprintf('SYLVESTER EQUATIONS : %s \n', num2str(nEquations))
fprintf('RANK REVEALING METRIC : %s \n', rank_revealing_metric)


% Get example polynomials
[fxy_exact, gxy_exact, hxy_exact,...
    uxy_exact,vxy_exact, wxy_exact,...
    dxy_exact,...
    m, m1, m2, n, n1, n2, o, o1, o2,...
    t_exact, t1_exact, t2_exact] = Examples_GCD_3Polys(ex_num);


fprintf('\n')
fprintf('----------------------------------------------------------------\n')
fprintf('Input Polynomials Degrees:\n')
fprintf('m  : %i \n', m)
fprintf('m1 : %i \n', m1)
fprintf('m2 : %i \n\n', m2)

fprintf('n  : %i \n', n)
fprintf('n1 : %i \n', n1)
fprintf('n2 : %i \n\n', n2)

fprintf('o  : %i \n', o)
fprintf('o1 : %i \n', o1)
fprintf('o2 : %i \n\n', o2)

fprintf('t  : %i \n', t_exact)
fprintf('t1 : %i \n', t1_exact)
fprintf('t2 : %i \n', t2_exact)
fprintf('----------------------------------------------------------------\n')
fprintf('\n')
fprintf('----------------------------------------------------------------\n')


% %
% %
% Add Noise

% Add noise to the coefficients of f and g
[fxy, ~] = AddVariableNoiseToPoly(fxy_exact, emin, emax);
[gxy, ~] = AddVariableNoiseToPoly(gxy_exact, emin, emax);
[hxy, ~] = AddVariableNoiseToPoly(hxy_exact, emin, emax);



% Get GCD d(x,y) and quotient polynomials u(x,y) and v(x,y)
lowerLimit = 1;
upperLimit = min([m,n,o]);
t_limits = [lowerLimit, upperLimit];

rank_range = [0 , 0];

% Get the GCD by my method
[fxy, gxy, hxy, dxy, uxy, vxy, wxy, t] = ...
    o_gcd_mymethod_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, t_limits, rank_range);



% Get Error in u(x,y), v(x,y) and d(x,y)
my_error.uxy = GetError(uxy, uxy_exact);
my_error.vxy = GetError(vxy, vxy_exact);
my_error.wxy = GetError(wxy, wxy_exact);
my_error.dxy = GetError(dxy, dxy_exact);

my_error.mean = geomean([my_error.uxy, my_error.vxy, my_error.wxy my_error.dxy]);


% Print the error in u(x,y), v(x,y) and d(x,y)
LineBreakLarge()
fprintf([mfilename ' : ' sprintf('Distance u(x,y) : %e \n', my_error.uxy)]);
fprintf([mfilename ' : ' sprintf('Distance v(x,y) : %e \n', my_error.vxy)]);
fprintf([mfilename ' : ' sprintf('Distance w(x,y) : %e \n', my_error.wxy)]);
fprintf([mfilename ' : ' sprintf('Distance d(x,y) : %e \n', my_error.dxy)]);
fprintf([mfilename ' : ' sprintf('Average : %e \n', my_error.mean)])
LineBreakLarge()

PrintToResultsFile(m, n, o, t, my_error)

end

function [dist] = GetDistance(fxy, gxy)
% Get Distance between two matrices
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)


%fxy = fxy./fxy(1,1);
%gxy = gxy./gxy(1,1);

fxy = fxy ./ norm(fxy);
gxy = gxy ./ norm(gxy);


sum_f = sum(sum(fxy));
sum_g = sum(sum(gxy));

if sum_f * sum_g < 0
   gxy = gxy.* -1; 
end

%fxy = MakePositive(fxy);
%gxy = MakePositive(gxy);

try
    dist = norm(fxy - gxy) ./ norm(fxy);
catch
    dist = 1000;
end

end


function fxy = MakePositive(fxy)

if fxy(1) < 0
    fxy = fxy.*-1;
end

end

function []= PrintToResultsFile(m,n,o,t,my_error)
%
% % Inputs
%
% m : (Int) Degree of polynomial f(x,y)
%
% n : (Int) Degree of polynomial g(x,y)
%
% o : (Int) Degree of polynomial h(x,y)
%
% t : (Int) Degree of GCD d(x,y)
%
% my_error : error.uxy, error.vxy, error.wxy and error.dxy

global SETTINGS

fullFileName = sprintf('Results/Results_o_gcd_3Polys.dat');

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
        
        % 19 FIELDS
        fprintf(fileID,'%s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s, %s \n',...
            datetime('now'),...
            SETTINGS.EX_NUM,...
            int2str(m),...
            int2str(n),...
            int2str(o),...
            int2str(t),...
            my_error.uxy,...
            my_error.vxy,...
            my_error.wxy,...
            my_error.dxy,...
            SETTINGS.MEAN_METHOD,...
            num2str(SETTINGS.BOOL_ALPHA_THETA),...
            SETTINGS.LOW_RANK_APPROX_METHOD,...
            num2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
            SETTINGS.APF_METHOD,...
            num2str(SETTINGS.APF_REQ_ITE),...
            SETTINGS.EMIN,...
            SETTINGS.EMAX,...
            SETTINGS.SYLVESTER_MATRIX_VARIANT...
            );
        % 19 arguments
        
    end

    function WriteHeader()
        
        headerString = ['DATE, EX_NUM, m, n, o, t, ERROR_UXY,' ...
            'ERROR_VXY, ERROR_WXY, ERROR_DXY, MEAN_METHOD,' ...
            'BOOL_ALPHA_THETA, LOW_RANK_APPROX_METHOD,'...
            'LRA_ITE, APF_METHOD, APF_ITE, error_min,' ...
            'error_max, Sylvester_Matrix_Variant \n'];
        
        fprintf(fileID,headerString);
    end

end





function [dist_fxy] = GetError(fxy,fxy_exact)
% GetError : Get distance between f(x,y) and exact form.
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of f(x,y) as computed.
%
% fxy_exact : (Matrix) Coefficients of f(x,y) exact.
%
% % Outputs.
%
% dist_fxy : Distance between f(x,y) and exact f(x,y)

dist_fxy = GetDistance(fxy_exact, fxy);

end


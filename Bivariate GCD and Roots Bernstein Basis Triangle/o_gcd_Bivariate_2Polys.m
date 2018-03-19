function [] = o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_type, rank_revealing_metric)
% O_GCD_BIVARIATE_2POLYS : Compute the GCD of two polynomials f(x,y) and g(x,y) where f(x,y)
% and g(x,y) are given in the Bernstein form, and are taken from the
% example file.
%
% % Inputs.
%
% ex_num : (String) Example Number
%
% el : (Float) Noise level - lower level
%
% em : (Float) Noise level - maximum level
%
% mean_method : (String) Options for mean calculation method
%               *Geometric Mean Matlab Method
%               *Geometric Mean My Method
%               *None
%
% bool_alpha_theta : (Boolean) Determines whether optimal values of alpha
% and theta are computed.
%           * true :  
%           * false :
%
% low_rank_approx_method : Options are (String)
%               * Standard STLN
%               * Standard SNTLN
%               * None
%
% apf_method : (String)
%       'None'
%       'Standard Linear APF'
%       'Standard Nonlinear APF'
%
% sylvester_matrix_type : (String)
%   * T
%   * DT
%   * TQ
%   * DTQ
%   * DTQ Denominator Removed
%
% rank_revealing_metric : (String)
%   * R1 Row Norms
%   * R1 Row Diagonals
%   * Minimum Singular Values
%   * Residuals
%
%
% % Examples
%
% >> o_gcd_Bivariate_2Polys('1', 1e-10, 1e-12, 'None', false, 'None', 'None', 'DTQ', 'Minimum Singular Values')
% >> o_gcd_Bivariate_2Polys('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'None', 'None', 'DTQ', 'Minimum Singular Values')
% >> o_gcd_Bivariate_2Polys('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'Standard STLN', 'None', 'DTQ', 'Minimum Singular Values')


% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    
    temp = emin;
    emin = emax;
    emax = temp;
    
end

% % Set global variables
SetGlobalVariables_GCD_2Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_type, rank_revealing_metric);


% Restore default path and add subfolders of current directory
restoredefaultpath
addpath(genpath(pwd))

% Print parameters to console
fprintf('INPUTS. \n')
fprintf('EXAMPLE NUMBER %s \n',ex_num)
fprintf('EMIN : %s \n',emin)
fprintf('EMAX : %s \n',emax)
fprintf('MEAN METHOD : %s \n', mean_method)
fprintf('PREPROCESSING : %s \n',num2str(bool_alpha_theta))
fprintf('RANK REVEALING METRIC : %s \n', rank_revealing_metric)
fprintf('LOW RANK METHOD : %s \n',low_rank_approx_method)
fprintf('APF METHOD : %s \n', apf_method)
fprintf('SYLVESTER FORMAT : %s \n', sylvester_matrix_type)


% %
% Get two polynomials f(x,y) and g(x,y) from example file.
[fxy, gxy, uxy_exact, vxy_exact, dxy_exact, m, n, t_exact] = ...
    Examples_GCD(ex_num);


% %
% Add noise to the coefficients

fxy = AddVariableNoiseToPoly(fxy, emin, emax);
gxy = AddVariableNoiseToPoly(gxy, emin, emax);

% %
% Get GCD by my method
limits_t = [0 min(m,n)];
rank_range = [0 0];

[fxy, gxy, dxy, uxy, vxy, t, rank_range] = ...
    o_gcd_mymethod_Bivariate_2Polys(fxy, gxy, m, n, limits_t, rank_range);


% Get Error in u(x,y), v(x,y) and d(x,y)
my_error.uxy = GetMatrixDistance(uxy_exact, uxy);
my_error.vxy = GetMatrixDistance(vxy_exact, vxy);
my_error.dxy = GetMatrixDistance(dxy_exact, dxy );

%average_error = mean([log10(my_error.uxy) log10(my_error.vxy) log10(my_error.dxy) ]);
average_error = geomean([(my_error.uxy) (my_error.vxy) (my_error.dxy) ]);
% Print the error in u(x,y), v(x,y) and d(x,y)
LineBreakLarge()
fprintf([mfilename ' : ' sprintf('Distance u(x,y) : %e \n', my_error.uxy)]);
fprintf([mfilename ' : ' sprintf('Distance v(x,y) : %e \n', my_error.vxy)]);
fprintf([mfilename ' : ' sprintf('Distance d(x,y) : %e \n', my_error.dxy)]);
fprintf([mfilename ' : ' sprintf('Average: %e \n', average_error)]);
LineBreakLarge()

PrintToResultsFile(m,n,t,my_error)

end



function [] = PrintToResultsFile(m, n, t, my_error)

global SETTINGS


fullFileName = sprintf('Results/Results_o_gcd.dat');

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
        
        % 15 FIELDS
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
            datetime('now'),...
            SETTINGS.EX_NUM,...
            int2str(m),...
            int2str(n),...
            int2str(t),...
            my_error.uxy,...
            my_error.vxy,...
            my_error.dxy,...
            SETTINGS.MEAN_METHOD,...
            num2str(SETTINGS.BOOL_ALPHA_THETA),...
            SETTINGS.LOW_RANK_APPROX_METHOD,...
            num2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
            SETTINGS.APF_METHOD,...
            num2str(SETTINGS.APF_REQ_ITE),...
            SETTINGS.EMIN,...
            SETTINGS.EMAX,...
            SETTINGS.SYLVESTER_MATRIX_VARIANT,...
            SETTINGS.RANK_REVEALING_METRIC...
            );
        
    end

    function WriteHeader()
        fprintf(fileID,'DATE, EX_NUM, m, n, t, ERROR_UXY, ERROR_VXY, ERROR_DXY, MEAN_METHOD, BOOL_ALPHA_THETA, LOW_RANK_APPROX_METHOD, LRA_ITE, APF_METHOD, APF_ITE, error_min, error_max, SUBRESULTANT_FORMAT, RANK_REVEALING_METRIC \n');
    end

end


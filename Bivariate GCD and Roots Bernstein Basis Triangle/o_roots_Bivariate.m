function [] = o_roots_Bivariate(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_variant, rank_revealing_metric, deconvolution_method_hxy, ...
    deconvolution_method_wxy, nEquations)
% O_ROOTS_BIVARIATE(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method, apf_method, sylvester_matrix_type, rank_revealing_method, deconvolution_method_hxy, deconvolution_method_wxy)
%
% Given an example number and set of parameters, obtain the roots of the
% example polynomial, where the polynomial is in the power basis
%
% % Inputs
%
% ex_num : (String) Example Number
%
% emin : (Float) Lower noise level
%
% emax : (Float) Upper noise level
%
% mean_method : (String)
%   * 'Geometric Mean Matlab Method'
%   * 'None'
%
% bool_alpha_theta : (Boolean)
%   * true :  Include Preprocessing
%   * false :  Exclude Preprocessing
%
% low_rank_approx_method : (String)
%   * 'Standard SNTLN' : Include SNTLN
%   * 'None'           : Exclude SNTLN
%
% sylvester_build_variant : (String) Variant of Sylvester matrix used
%   * 'T'
%   * 'DT'
%   * 'TQ'
%   * 'DTQ'
%   * 'DTQ Denominator Removed'
%    
%
%
%
%   
%
% % Examples
% >> o_roots_Bivariate('1', 1e-12, 1e-10, 'None', false, 'None', 'None', 'DTQ', 'Minimum Singular Values', 'Batch', 'Batch')
% >> o_roots_Bivariate('1', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'Standard STLN', 'Standard APF', 'DTQ', 'Minimum Singular Values', 'Batch', 'Batch')



% Restore defaults and add subfolders
restoredefaultpath
addpath(genpath(pwd));


SetGlobalVariables_Roots(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_variant, rank_revealing_metric, ...
    deconvolution_method_hxy, deconvolution_method_wxy, nEquations);


% %
% %
% Get Example

% Given the example number get the polynomial coefficients of f(x,y) and
% the total degree 'm' of f(x,y).
[fxy_exact, arr_fxy_exact, arr_hxy_exact, arr_wxy_exact, m] = Examples_Roots(ex_num);

% Add noise to the coefficients of polynomial f(x,y)
[fxy, ~] = AddNoiseToPoly(fxy_exact, emin);


% Get the array of polynomials f_{i}(x,y), h_{i}(x,y), w_{i}(x,y)        
[arr_fxy, arr_hxy, arr_wxy, ~] = o_roots_mymethod_newmethod(fxy, m);
        

 

% Compute the distance between the array or exact polynomials f_{i}(x,y)
% and the computed polynomials f_{i}(x,y)
vErrors_arr_fxy = CompareArrays(arr_fxy_exact, arr_fxy);

% Compute the distance between the array of exact polynomials h_{i}(x,y)
% and the computed polynomials h_{i}(x,y)
vErrors_arr_hxy = CompareArrays(arr_hxy_exact, arr_hxy);

% Compute the distance between the array of exact polynomials w_{i}(x,y)
% and the computed polynomials w_{i}(x,y)
vErrors_arr_wxy = CompareArrays(arr_wxy_exact, arr_wxy);


% Plot Errors in the set of polynomials {f_{i}(x,y)}, 
% the set of polynomials {h_{i}(x,y)} and 
% the set of polynomials {w_{i}(x,y)} 
PlotErrors({vErrors_arr_fxy, vErrors_arr_hxy, vErrors_arr_wxy},...
    {'$\epsilon f_{i}(x,y)$', ...
    '$\epsilon h_{i}(x,y)$', ...
    '$\epsilon w_{i}(x,y)$'});


my_error.arr_fxy = mean(vErrors_arr_fxy);
my_error.arr_hxy = mean(vErrors_arr_hxy);
my_error.arr_wxy = mean(vErrors_arr_wxy);

% Print Errors
LineBreakLarge()
fprintf('Average Error f_{i}(x,y) : %e \n', my_error.arr_fxy);
fprintf('Average Error h_{i}(x,y) : %e \n', my_error.arr_hxy);
fprintf('Average Error w_{i}(x,y) : %e \n', my_error.arr_wxy);
LineBreakLarge()


PrintToResultsFile(my_error)


end

function [vDistance] = CompareArrays(arr_fxy_exact, arr_fxy)


% Get length of array
nPolysArray = length(arr_fxy_exact);

% Store distance
vDistance = zeros(nPolysArray,1);

for i = 1 : 1 : nPolysArray
    
    fxy_exact = arr_fxy_exact{i};
    fxy_comp = arr_fxy{i};

    vDistance(i) = GetMatrixDistance(fxy_exact, fxy_comp);
end



end


function PlotErrors(arrErrors, arrStrings)
%
% % Input
%
% arrErrors : (Array of Vectors) Each cell of the array contains a vector
% of errors
%
%
% arrStrings : (Array of Strings)
figure_name = 'Errors';
figure('Name',figure_name)

hold on
for i = 1 : 1 : 3
   
    
    vErrors = arrErrors{i};
    nErrors = length(vErrors);
    x_vec = 1 : 1 : nErrors;
    strName = arrStrings{i};
    plot(x_vec, log10(vErrors), '-o', 'DisplayName', strName, 'LineWidth', 2);
    
end

l = legend(gca,'show');
set(l,'Interpreter', 'latex');

grid on
box on

hold off

end


function [] = PrintToResultsFile(my_error)

global SETTINGS

% Get datetime and filename
%v = datevec(now);
%fullFileName = sprintf('Results/Results_o_gcd_%s-%s-%s.txt',num2str(v(1)), num2str(v(2)), num2str(v(3)));

fullFileName = sprintf('Results/Results_o_roots.dat');

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
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
            datetime('now'),...
            SETTINGS.EX_NUM,...
            my_error.arr_fxy,...
            my_error.arr_hxy,...
            my_error.arr_wxy,...
            SETTINGS.MEAN_METHOD,...
            num2str(SETTINGS.BOOL_ALPHA_THETA),...
            SETTINGS.LOW_RANK_APPROX_METHOD,...
            num2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
            SETTINGS.APF_METHOD,...
            num2str(SETTINGS.APF_REQ_ITE),...
            SETTINGS.EMIN,...
            SETTINGS.EMAX,...
            SETTINGS.SYLVESTER_MATRIX_VARIANT,...
            SETTINGS.RANK_REVEALING_METRIC,...
            SETTINGS.DECONVOLUTION_METHOD_HXY,...
            SETTINGS.DECONVOLUTION_METHOD_WXY...
            );
        
    end

    function WriteHeader()
        myString = ['DATE, EX_NUM, ERROR_FXY, ERROR_HXY, ERROR_WXY,' ...
            'MEAN_METHOD, BOOL_ALPHA_THETA, LOW_RANK_APPROX_METHOD, LRA_ITE, '...
            'APF_METHOD, APF_ITE, error_min, error_max, SUBRESULTANT_FORMAT, '...
            'RANK_REVEALING_METRIC, DECONVOLUTION_METHOD_HXY, DECONVOLUTION_METHOD_WXY \n'];
        
        fprintf(fileID, myString);
    end

end

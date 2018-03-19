
function [] = PrintRootsToFile(arr_RootFindingMethod, arr_BackwardErrors, ...
    arr_ForwardErrors, vError_arr_fx, vError_arr_hx, vError_arr_wx)
% Print results of root finding computation to a text file
%
% % Inputs
%
% arr_RootFindingMethod : (Array of Strings) Array of the names of the root
% finding methods used
%
% arr_ForwardErrors : 
%
% arr_BackwardErrors : 
%
% arr_errors : (Array of Floats)
%
% vErrors_arr_fx : (Vector) Vector of errors where the ith entry
% corresponds to the error in f_{i}(x)
%
% vErrors_arr_hx : (Vector) Vector of errors where the ith entry
% corresponds to the error in h_{i}(x)
%
% vErrors_arr_wx : (Vector) Vector of errors where the ith entry
% corresponds to the error in w_{i}(x)


global SETTINGS


nMethods = length(arr_RootFindingMethod);

fullFileName = sprintf('Results/Results_o_roots.dat');

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    
    for i = 1 : 1 : nMethods
        method_name = arr_RootFindingMethod{i};
        method_BackwardError = arr_BackwardErrors{i};
        method_ForwardError = arr_ForwardErrors{i};
        
        error_fx = norm(vError_arr_fx);
        error_hx = norm(vError_arr_hx);
        error_wx = norm(vError_arr_wx);
        
        % Only print results if my method
        if (strcmp(method_name, 'My Method'))
            WriteNewLine(method_name, method_BackwardError, ...
                method_ForwardError, error_fx, error_hx, error_wx);
        end
        
        
    end
    
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    for i = 1 : 1 : nMethods
        
        method_name = arr_RootFindingMethod{i};
        method_BackwardError = arr_BackwardErrors{i};
        method_ForwardError = arr_ForwardErrors{i};
        
        % Only print results if my method
        if (strcmp(method_name, 'My Method'))
            WriteNewLine(method_name, method_BackwardError, method_ForwardError)
        end
        
    end
    fclose(fileID);
    
end



    function WriteNewLine(method_name, myBackwardError, myForwardError, error_fx, error_hx, error_wx)
        % Write a new line of the text file
        
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
            datetime('now'),...
            SETTINGS.EX_NUM,...
            SETTINGS.MEAN_METHOD,...
            num2str(SETTINGS.BOOL_ALPHA_THETA),...
            SETTINGS.EMIN,...
            SETTINGS.EMAX,...
            SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
            num2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
            SETTINGS.APF_METHOD,...
            num2str(SETTINGS.APF_REQ_ITE),...
            num2str(SETTINGS.BOOL_LOG),...
            SETTINGS.SYLVESTER_MATRIX_VARIANT,...
            SETTINGS.GCD_COEFFICIENT_METHOD,...
            method_name,...
            SETTINGS.DECONVOLUTION_METHOD_HX,...
            SETTINGS.DECONVOLUTION_METHOD_WX,...
            num2str(SETTINGS.PREPROC_DECONVOLUTIONS),...
            myBackwardError,...
            myForwardError,...
            error_fx,...
            error_hx,...
            error_wx...
            );
    end


    function WriteHeader()
        % If the file doesnt already exist, write a header to the text file
        strHeaders = ['DATE, EX_NUM, MEAN_METHOD, BOOL_ALPHA_THETA, EMIN,'...
            'EMAX, LOW_RANK_APPROX_METHOD, LOW_RANK_ITE, APF_METHOD, '...
            'APF_ITE, BOOL_LOG, SYLVESTER_MATRIX_VARIANT, GCD_METHOD,'...
            'METHOD_NAME, DECONVOLUTION_METHOD HX, '...
            'DECONVOLUTION_METHOD WX, DECONVOLUTION_PREPROC,' ...
            'FORWARD_ERROR, BACKWARD_ERROR, ERROR_FX, ERROR_HX, ERROR_WX \n'];
        
        fprintf(fileID, strHeaders);
    end


end

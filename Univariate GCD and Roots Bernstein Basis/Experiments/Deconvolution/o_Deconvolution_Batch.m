function [] = o_Deconvolution_Batch
%
% This function performs a batch of deconvolution problems for various
% noise levels, both with and without preprocessing


% Initialise an array of example numbers
arr_ex_num = {'1', '2', '3', '4', '5', '6','7','8'};

% Initialise an array of noise levels
arr_noise = {1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-12};

% Initialise an array to store true and false for preprocessing option of
% the deconvolution method
arr_bool_preproc = {true, false};


% Set the filename for the log file which will log any errors or successes
% of the deconvolution method
log_filename = 'log-deconvolution.txt';


parfor i1 = 1:1:length(arr_ex_num)
    
    % Get example number
    ex_num = arr_ex_num{i1};
    
    for i2 = 1:1:length(arr_noise)
    
        % Get minimum noise level
        emin = arr_noise{i2};
        
        for i3 = 1:1:length(arr_bool_preproc)
            
            % Get bool_preproc from the preproc array. This determines 
            % whether preprocessing is used in the deconvolution function
            bool_preproc = arr_bool_preproc{i3};
            
            
            try
                
                close all;
                clc;
                o_Deconvolution(ex_num, emin, bool_preproc)
                
                % Print success to log file
                fileId = fopen(log_filename, 'a');
                
                fprintf(fileId, '%s %s \n',...
                    datetime('now'),...
                    'Success'...
                );
                
                fclose(fileId);
                
            catch err
                % Print failure to log file
                fileId = fopen(log_filename, 'a');
                
               
                
                fprintf(fileId, '%s %s \n', ...
                    datetime('now'),...
                    getReport(err)...
                );
                
                fclose(fileId);
           end
        end
    end
end

    


end
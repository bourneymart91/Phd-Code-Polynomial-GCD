function [] = o_Deconvolution_batch
% Perform a batch of deconvolutions
%

ex_num_arr = {'1','2','3','4'};

parfor i1 = 1:1:length(ex_num_arr)
    ex_num = ex_num_arr{i1};
    
    
    try
        close all;
        clc;
        o_Deconvolution(ex_num);
        
        % Print success to log file
        fileId = fopen('log.txt','a');
        fprintf(fileId,'%s','success \n');
        fclose(fileId);
        
    catch err
        % Print failure to log file
        fileId = fopen('log.txt','a');
        fprintf(fileId,'%s \n\n\n',getReport(err));
        fclose(fileId);
    end
    
end

end
function [] = o_deconvolution_Batch()
% Perform a batch of deconvolution examples
% 
% % Inputs


global SETTINGS

SETTINGS.PLOT_GRAPHS_DECONVOLUTION = false;

% Get Example Numbers
arr_ex_num = {'1', '2', '3'};

% Get noise levels
arr_noise = {1e-6, 1e-8, 1e-10, 1e-12};


arr_bool_preproc = {true, false};

% For each example number
for i1 = 1 : 1 : length(arr_ex_num)
    ex_num = arr_ex_num{i1};

    % For each noise level
    for i2 = 1 : 1 : length(arr_noise)
        emin = arr_noise{i2};
    
        % With and without preprocessing
        for i3 = 1 : 1 : length(arr_bool_preproc)
            bool_preproc = arr_bool_preproc{i3};
            
            o_deconvolution(ex_num, emin, bool_preproc)
            
        end
    end
end


end
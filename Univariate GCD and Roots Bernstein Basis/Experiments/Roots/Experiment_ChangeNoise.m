function [] = Experiment_ChangeNoise(ex_num)
% This experiment repeats a root finding example for a variety of noise
% levels
%
% % Inputs
%
% ex_num : (String) Example number

close all;
clc;

% Set minimum noise level
emin = 1e-12;

% Get array of maximum noise levels
arrEmax = {1e-10, 1e-8, 1e-6};

% Set preprocessing related variables
switch bool_preproc
    case true
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
    case false
        mean_method = 'None';
        bool_alpha_theta = false;
end


low_rank_approx_method = 'None';
apf_method = 'None';

% Set the sylvester matrix variant to be used.
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_matrix_variant = 'DTQ';

% Method used to determine the degree of the GCD
% 'R1 Row Norms',
% 'R1 Row Diagonals',
% 'Minimum Singular Values',
% 'Normalised Minimum Singular Values'
rank_revealing_metric = 'Minimum Singular Values';


arrDeconvolution_method = {'Separate', ...
    'Batch', ...
    'Batch With STLN', ...
    'Batch Constrained', ...
    'Batch Constrained With STLN'};

arrDeconvolution_preproc = {true, false};

arrDeconvolution_method_wx = {'Batch', ...
    'Separate'};


parfor i = 1:1:length(arrDeconvolution_method)
    for i2 = 1:1: length(arrDeconvolution_preproc)
        for i3 = 1: 1: length(arrDeconvolution_method_wx)
            for i4 = 1 : 1 : length(arrEmax)
                
                emax = arrEmax{i4};
                deconvolution_method_hx = arrDeconvolution_method{i};
                deconvolution_preproc = arrDeconvolution_preproc{i2};
                deconvolution_method_wx = arrDeconvolution_method_wx{i3};
                
                o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
                    low_rank_approx_method, apf_method, sylvester_matrix_variant, ...
                    rank_revealing_metric, deconvolution_method_hx, ...
                    deconvolution_method_wx, deconvolution_preproc)
                
            end
        end
    end
end
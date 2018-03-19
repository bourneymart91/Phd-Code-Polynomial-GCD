function [] = Experiment1(ex_num)
%
% % Input
%
% ex_num : (String) Example Number
%


close all;
clc;

% Set upper and lower noise level
emin = 1e-10;
emax = 1e-12;

% Set preprocessing related variables
bool_preproc = true;

switch bool_preproc
    case true
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
    case false
        mean_method = 'None';
        bool_alpha_theta = false;        
end

% Set other variables
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

%deconvolution_method_hx = 'Separate';
deconvolution_method_wx = 'Batch';

% For each deconvolution method
for i1 = 1 : 1 : length(arrDeconvolution_method)
    
    for i2 = 1:1: length(arrDeconvolution_preproc)
        
        % Set method for deconvolution of polynomials h_{i}(x)
        deconvolution_method_hx = arrDeconvolution_method{i1};
        
        % Set boolean value true or false as to whether polynomials
        % h_{i}(x) are preprocessed
        deconvolution_preproc = arrDeconvolution_preproc{i2};
        
        % Compute roots
        o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
            low_rank_approx_method, apf_method, sylvester_matrix_variant, ...
            rank_revealing_metric, deconvolution_method_hx, ...
            deconvolution_method_wx, deconvolution_preproc)
        
        
        
    end
end




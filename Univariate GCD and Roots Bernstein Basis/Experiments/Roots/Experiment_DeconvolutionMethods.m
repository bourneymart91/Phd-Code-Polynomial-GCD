function [] = Experiment_DeconvolutionMethods(ex_num)
% This experiment considers the various deconvolution methods when
% deconvolving the set of polynomials f_{i}(x) to obtain the set of
% polynomials h_{i}(x)
%
% % Inputs
%
% ex_num : (String) Example number
%
%
% % Examples
%
% >> Experiment_DeconvolutionMethods('1')


close all;
clc;


%
% -------------------------------------------------------------------------
% Constants

% Set noise level
emin = 1e-10;
emax = 1e-10;

% Set variables relating to preprocessing
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

% % Sylvester Build Method
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_matrix_variant = 'DTQ';

%apf_method
%   'None'
apf_method = 'None';

% Method used to determine the degree of the GCD
% 'R1 Row Norms',
% 'R1 Row Diagonals',
% 'Minimum Singular Values',
% 'Normalised Minimum Singular Values'
rank_revealing_metric = 'Minimum Singular Values';


% % Low Rank Approximation Method
% 'None'
% 'Standard STLN';
% 'Standard SNTLN'
low_rank_approx_method = 'None';



% Deconvolution method
% 'Separate'
% 'Batch'
deconvolution_method_wx = 'Batch';

%
% -------------------------------------------------------------------------
% Variables


%
deconvolution_method_hx_arr = {...
    'Separate', ...
    'Batch', ...
    'Batch with STLN', ...
    'Batch Constrained', ...
    'Batch Constrained with STLN'};



% -------------------------------------------------------------------------



% 
bool_deconvolution_preproc = false;

for i = 1 : 1 : length(deconvolution_method_hx_arr)
    
    deconvolution_method_hx = deconvolution_method_hx_arr{i};
    
    
    o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
        low_rank_approx_method, apf_method, sylvester_matrix_variant, ...
        rank_revealing_metric, deconvolution_method_hx, ...
        deconvolution_method_wx, bool_deconvolution_preproc)
    
end



% -------------------------------------------------------------------------

% Set preprocessing in the deconvolution problem = true
bool_deconvolution_preproc = true;


deconvolution_method_hx_arr = {...
    'Separate', ...
    'Batch', ...
    'Batch with STLN', ...
    'Batch Constrained', ...
    'Batch Constrained with STLN'};


parfor i = 1 : 1 : length(deconvolution_method_hx_arr)
    
    deconvolution_method_hx = deconvolution_method_hx_arr{i};
    
    
    o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
        low_rank_approx_method, apf_method, sylvester_matrix_variant, ...
        rank_revealing_metric, deconvolution_method_hx, ...
        deconvolution_method_wx, bool_deconvolution_preproc)
    
end


sameaxes()

end



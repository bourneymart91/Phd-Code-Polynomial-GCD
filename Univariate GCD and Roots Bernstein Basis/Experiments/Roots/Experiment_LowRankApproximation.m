function [] = Experiment_LowRankApproximation(ex_num, bool_preproc)
%
%
%
% >> Experiment_LowRankApproximation('1')
% >> Experiment_LowRankApproximation('15')
close all; clc;


%
% -------------------------------------------------------------------------
% Constants

emin = 1e-10;
emax = 1e-8;

% Set preprocessing related variables
switch bool_preproc
    case true
        
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case false
        
        mean_method = 'None';
        bool_alpha_theta = false;
        
end

% apf_method
%   'None'

apf_method = 'None';

% Set the sylvester matrix variant to be used.
%   'T'
%   'DT'
%   'TQ'
%   'DTQ'
sylvester_matrix_variant = 'DTQ';

% Method used to determine the degree of the GCD
%   'R1 Row Norms',
%   'R1 Row Diagonals',
%   'Minimum Singular Values',
%   'Normalised Minimum Singular Values'
rank_revealing_metric = 'Minimum Singular Values';

% deconvolution_method_hx
%   'Separate'
%   'Batch'
%   'Batch with STLN'
%   'Batch Constrained'
%   'Batch Constrained with STLN'
deconvolution_method_hx = 'Batch Constrained';

% deconvolution_method_wx
%   'Separate'
%   'Batch'
deconvolution_method_wx = 'Batch';

% bool_deconvolution_preproc
%   true  :
%   false :
bool_deconvolution_preproc = true;


arrLowRankApproximationMethods = {'None', 'Standard STLN', 'Standard SNTLN'};


for i = 1 : 1 : length(arrLowRankApproximationMethods
    
    
    low_rank_approx_method = arrLowRankApproximationMethods{i};
    
    
    
    o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
        low_rank_approx_method, apf_method, sylvester_matrix_variant, ...
        rank_revealing_metric, deconvolution_method_hx, ...
        deconvolution_method_wx, bool_deconvolution_preproc)
end

end



function [] = Experiment0(ex_num, bool_preproc)
% Given an example number and a boolean value, compute the roots of the
% univariate polynomial f(x)
%
%
% >> Experiment0('1')
%
% % Inputs
%
% ex_num : (String) Example Number
%
% bool_preproc : (Boolean) True or false determines whether polynomials are
% preprocessed.
%
%
% >> Experiment0('1', true)

close all; 
clc;

% Set upper and lower bound for noise level
emin = 1e-4;
emax = 1e-4;

% Set preprocessing related variables
switch bool_preproc
    case true
        
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case false
        
        mean_method = 'None';
        bool_alpha_theta = false;
        
end


% Set other variables
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

% % Low Rank Approximation Method
% 'None'
% 'Standard STLN'
% 'Standard SNTLN'
low_rank_approx_method = 'None';


% % Deconvolution Method
% 'Separate'
% 'Batch'
% 'Batch with STLN'
% 'Batch Constrained'
% 'Batch Constrained with STLN'
deconvolution_method_hx = 'Batch Constrained with STLN';

% % Deconvolution Method
% 'Batch'
% 'Separate'
deconvolution_method_wx = 'Batch';


% Determine whether to preprocess the set of polynomials in the
% deconvolution problem
deconvolution_preproc = true;

% Compute the roots of f(x)
o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_matrix_variant, ...
    rank_revealing_metric, deconvolution_method_hx, ...
    deconvolution_method_wx, deconvolution_preproc)


end


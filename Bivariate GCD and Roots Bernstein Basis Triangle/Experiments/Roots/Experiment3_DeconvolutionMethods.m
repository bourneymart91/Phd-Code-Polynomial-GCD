function [] = Experiment3_DeconvolutionMethods(ex_num, bool_preproc)
% Experiment with various deconvolution methods
%
% % Inputs
%
% ex_num : (String) Example number
%
% bool_preproc : (Boolean)
%
%


close all; clc;


emin = 1e-10;
emax = 1e-8;


switch bool_preproc
    
    case true
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case false
        mean_method = 'None';
        bool_alpha_theta = false;
end

% Method used to add structure to t-th subresultant matrix such that a low
% rank approximation is obtained.
%   'None'
%   'Standard STLN'
low_rank_approx_method = 'None';

% Method used to add structure to t-th [C(u) ; C(v)] matrix
apf_method = 'None';

% Set variant of subresultant matrix
%   'T'
%   'DT'
%   'TQ'
%   'DTQ'
sylvester_matrix_variant = 'DTQ';

% nEquations determines structure of the Sylvester matrix 
%   '2' : 2 x 3 partitioned structure
%   '3' : 3 x 3 partitioned structure
nEquations = '2';

% Choose method for determining the degree of the GCD
%   'Minimum Singular Values'
%   ''
rank_revealing_method = 'Minimum Singular Values';




% Choose method for deconvolution of the set of polynomials f_{i}(x,y) to
% obtain the set of polynomials h_{i}(x,y)
%   'Separate'
%   'Batch'
%   'Batch with STLN'
%   'Batch Constrained'
%   'Batch Constrained with STLN'
arr_deconvolution_method_hxy = {'Separate', 'Batch',  'Batch With STLN', ...
    'Batch Constrained', 'Batch Constrained With STLN'};


% Chose method of deconvolution for the 
deconvolution_method_wxy = 'Batch';




nDeconvolutionMethods = length(arr_deconvolution_method_hxy);

for i = 1 : 1 : nDeconvolutionMethods
    
    deconvolution_method_hxy = arr_deconvolution_method_hxy{i};
    
    
% Compute factorisation
o_roots_Bivariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_matrix_variant, ...
    rank_revealing_method, deconvolution_method_hxy, deconvolution_method_wxy, nEquations)
end

end
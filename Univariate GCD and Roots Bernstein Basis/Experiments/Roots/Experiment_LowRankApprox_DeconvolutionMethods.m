function [] = Experiment_LowRankApprox_DeconvolutionMethods(ex_num, bool_preproc)
%
% % Inputs
%
% ex_num : (String) Example number
%
% bool_preproc : (Boolean) Determine whether preprocessing is used in the
% computation of the degree of the GCD in each GCD computation as part of
% the factorisation algorithm
%
% % Examples
%
% >> Experiment_LowRankApprox_DeconvolutionMethods('10')

close all;
clc;


%
% -------------------------------------------------------------------------
% Constants

% Upper and lower noise level
emin = 1e-4;
emax = 1e-4;

% Preprocessing related variables
switch bool_preproc
    
    case true
        
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case false
        
        mean_method = 'None';
        bool_alpha_theta = false;
        
end


% Other variables 

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



%
% -------------------------------------------------------------------------
% Variables


deconvolution_method_wx = 'Batch';
bool_deconvolution_preproc = true;

%arr_deconvolution_method_hx = {'Separate', 'Batch', 'Batch With STLN', 'Batch Constrained', 'Batch Constrained With STLN'};
arr_deconvolution_method_hx = {'Separate', ...
    'Batch'};

arr_low_rank_approx_method = {'None', ...
    'Standard SNTLN'};



for i1 = 1 : 1 : length(arr_deconvolution_method_hx)
    
    deconvolution_method_hx = arr_deconvolution_method_hx{i1};
    
    for i2 = 1 : 1 : length(arr_low_rank_approx_method)
        
        low_rank_approx_method = arr_low_rank_approx_method{i2};
        
        
        try
        [epsilon_fi(i1,i2), epsilon_hi(i1,i2), epsilon_wi(i1,i2)] = ...
            o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
            low_rank_approx_method, apf_method, sylvester_matrix_variant, ...
            rank_revealing_metric, deconvolution_method_hx, ...
            deconvolution_method_wx, bool_deconvolution_preproc)
        
        
        catch
        
        end
        
    end
    
end

display(epsilon_fi)

display(epsilon_hi)

display(epsilon_wi)

sameaxes()


end



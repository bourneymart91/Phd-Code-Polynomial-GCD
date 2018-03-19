function [] = Experiment2Preprocessing_3Polys(ex_num, ex_num_variant)
% >> Experiment2Preprocessing_3Polys(ex_num, ex_num_variant)
% 
% This experiment considers the computation of the GCD of three univariate
% polynomials by Sylvester subresultant matrix based methods, where the
% matrices may or may not be preprocessed. It is typically shown that
% preprocessing yields improved results.
%
% % Inputs
%
% ex_num : (String) Example number
%
% ex_num_variant : (String) 'a', 'b' or 'c'
%
%
% % Outputs
%
% Results are printed to the .dat file called "Results_o_gcd_3Polys.dat"
% contained in the "Results" folder.
%
% >> Experiment2Preprocessing_3Polys('1', 'a')


% Clear all
close all; clc;


% Constants ===============================================================

% Set minimum and maximum level of componentwise noise
emin = 1e-7;
emax = 1e-4;

% Constants

% Low Rank Approximation Method 
%   'None'
%   'Standard STLN'
%   'Standard SNTLN'
low_rank_approx_method = 'None';


apf_method = 'None';

% Sylvester Matrix Variant
% Set the sylvester matrix variant to be used. 
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_matrix_variant = 'DTQ';

% rank_revealing_metric
%   'Normalised Minimum Singular Values'
%   'Minimum Singular Values'
rank_revealing_metric = 'Normalised Minimum Singular Values';



global SETTINGS

SETTINGS.SCALING_METHOD = 'NONE';

% nEquations 
%   '2' : Subresutlant matrices have a 2 by 3 partitioned structure
%   '3' : Subresultant matrices have a 3 by 3 partitioend structure
nEquations = '2';



% Variables ===============================================================

% Exclude Preprocessing
mean_method = 'None';
bool_alpha_theta = false;



o_gcd_Univariate_3Polys(ex_num, ex_num_variant, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_variant, nEquations, rank_revealing_metric)





% Include Preprocessing

bool_alpha_theta = true;
mean_method = 'Geometric Mean Matlab Method';

% Must set a scaling method
arrScaling_method = {'lambda_rho'};

for i = 1 : 1 : length(arrScaling_method)
    
    SETTINGS.SCALING_METHOD = arrScaling_method{i};
    o_gcd_Univariate_3Polys(ex_num, ex_num_variant, emin, emax, mean_method, ...
        bool_alpha_theta, low_rank_approx_method, apf_method, ...
        sylvester_matrix_variant, nEquations, rank_revealing_metric)
    
end

end
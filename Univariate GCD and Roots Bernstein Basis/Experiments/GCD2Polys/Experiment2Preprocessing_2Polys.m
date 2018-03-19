function [] = Experiment2Preprocessing_2Polys(ex_num)
% Experiment performs GCD computation with and without preprocessing
%
%
% % Input
%
% ex_num : (String) Example number



% Constants : 
% Variable : Preprocessing 


close all; clc;

% Set experiment constants ----------------------------------------

% Set upper and lower noise limit
emin = 1e-6;
emax = 1e-8;


% low_rank_approx_method 
%   'None'
%   'Standard STLN'
%   'Standard SNTLN'
low_rank_approx_method = 'None';

% apf_method 
%   'None'
apf_method = 'None';

% Sylvester_matrix_variant 
%   'T'
%   'DT'
%   'TQ'
%   'DTQ'
sylvester_matrix_variant = 'DTQ';

% rank_revealing_metric
%   'Minimum Singular Values'
%   ''
rank_revealing_metric = 'Minimum Singular Values';

% Variables ------------------------------------------

% Perform GCD computation without preprocessing
mean_method = 'None';
bool_alpha_theta = false;

o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_matrix_variant, rank_revealing_metric) ;

% Perform GCD computation with preprocessing
mean_method = 'Geometric Mean Matlab Method';
bool_alpha_theta = true;

o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_matrix_variant, rank_revealing_metric) ;



end
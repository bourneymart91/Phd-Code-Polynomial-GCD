Experiment6RankRevealingMetric_2Polys(ex_num)
% This experiment considers the optimal method for the computation of the
% degree of the GCD by numerical rank of the set of subresultant matrices
%
% % Inputs
%
% ex_num : (String) Example Number


close all;
clc;

% Set upper and lower noise level
el = 1e-8;
eu = 1e-6;

% Determine whether preprocessing is used
bool_preprocessing = true;
switch bool_preprocessing
    case 1
        
        mean_method = 'Geometric Mean Matlab Method';
        bool_preproc = true;
    case 0
        mean_method = 'None';
        bool_preproc = false;
end

% Set the sylvester matrix variant to be used. 
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_matrix_variant = 'DTQ';

% 'None'
% 'Standard STLN'
% 'Standard SNTLN'
low_rank_approx_method = 'None';

%
apf_method = 'None';

% Method used to determine the degree of the GCD
% 'R1 Row Norms',
% 'R1 Row Diagonals',
% 'Minimum Singular Values',
% 'Normalised Minimum Singular Values'
rank_revealing_metric_arr = {'Minimum Singular Values', ...
    'Normalised Minimum Singular Values', ...
    'R1 Row Diagonals',...
    'R1 Row Norms'};

for i = 1 : 1 : length(rank_revealing_metric_arr)
    
    rank_revealing_metric = rank_revealing_metric_arr{i};
    
    o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, ...
        low_rank_approx_method, apf_method, sylvester_matrix_variant, ...
        rank_revealing_metric) ;
end
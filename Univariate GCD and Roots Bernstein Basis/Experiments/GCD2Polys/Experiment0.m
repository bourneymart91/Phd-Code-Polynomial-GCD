function [] = Experiment0(ex_num, bool_preproc)
% Experiment with optimal low rank approximation method for the t-th
% subresultant matrix. Where 't' is the degree of the GCD.
%
% % Inputs
%
% ex_num : (String) Example number
%
% bool_preproc : (Boolean) 



close all; clc;


% Set upper and lower noise level
el = 1e-12;
eu = 1e-10;


% Method used to determine the degree of the GCD
% 'R1 Row Norms',
% 'R1 Row Diagonals',
% 'Minimum Singular Values',
% 'Normalised Minimum Singular Values'
rank_revealing_metric = 'Minimum Singular Values';

% Set preprocessing related variables


if (bool_preproc == true)
    
    mean_method = 'Geometric Mean Matlab Method';
    bool_alpha_theta = true;
    
else
    mean_method = 'None';
    bool_alpha_theta = false;
    
end

% APF method
apf_method = 'None';

% Set the sylvester matrix variant to be used. 
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_matrix_variant = 'DTQ';


% Get low rank approximation methods
low_rank_approx_method = 'None';

    
o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_matrix_variant, rank_revealing_metric) ;




end


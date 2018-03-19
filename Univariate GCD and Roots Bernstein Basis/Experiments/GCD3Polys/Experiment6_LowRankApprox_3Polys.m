function [] = Experiment6_LowRankApprox_3Polys(ex_num, bool_preproc)
% Experiment considers low rank approximation methods
%
% % Inputs
%
%
% ex_num : (String) Example number
%
% bool_preproc : (Boolean) Preprocessing
%
%



close all; clc;

% Noise levels
emin = 1e-10;
emax = 1e-8;

% nEquations : Determines the structure of the three-polynomial subresulant
% matrix 
%       '2' : Subresultant matrices have a 2 by 3 partitioned structure
%       '3' : Subresutlant matrices have a 3 by 3 partitioned structure
nEquations = '2';

% apf_method 
%   'None'
apf_method = 'None';


% Set the sylvester matrix variant to be used.
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_matrix_variant = 'DTQ';


% rank_revealing_metric
%   'Minimum Singular Values'
%   'Normalised Minimum Singular Values'
rank_revealing_metric = 'Minimum Singular Values';



global SETTINGS

SETTINGS.SCALING_METHOD = 'lambda_mu_rho';

switch bool_preproc
    case 1
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case 0
        mean_method = 'None';
        bool_alpha_theta = false;
        
    otherwise
        error('err')
end


low_rank_approx_method = 'None';

o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta,...
    low_rank_approx_method, apf_method, sylvester_matrix_variant, nEquations, ...
    rank_revealing_metric)



low_rank_approx_method = 'Standard STLN';

o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta,...
    low_rank_approx_method, apf_method, sylvester_matrix_variant, nEquations, ...
    rank_revealing_metric)


end



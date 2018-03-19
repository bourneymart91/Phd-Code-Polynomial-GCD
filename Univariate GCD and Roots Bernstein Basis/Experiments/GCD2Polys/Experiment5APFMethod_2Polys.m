function [] = Experiment5APFMethod_2Polys(ex_num)
% Experiment with the optimal APF method. Where the APF method is defined
% as...
%
% Inputs
%
% ex_num : (String) Example number


close all;
clc;


% Set upper and lower noise level.
el = 0;
eu = 0;

switch bool_preprocessing
    case 1
        mean_method = 'Geometric Mean Matlab Method';
        bool_preproc = true;
        
    case 0
        mean_method = 'None';
        bool_preproc = false;
        
end


% Sylvester matrix variant 
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_matrix_variant = 'T';

% % Method used to determine the degree of the GCD
% R1 Row Norms,
% R1 Row Diagonals,
% Minimum Singular Values,
% Normalised Minimum Singular Values
rank_revealing_metric = 'Minimum Singular Values';

% Set method for computing the low rank approximation of the t-th
% subresultant matrix 
% None
% Standard STLN
% Standard SNTLN
low_rank_approx_method = 'None';


apf_method = 'Last Non-zero Row';

o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, ...
    low_rank_approx_method, apf_method, sylvester_matrix_variant, ...
    rank_revealing_metric) ;




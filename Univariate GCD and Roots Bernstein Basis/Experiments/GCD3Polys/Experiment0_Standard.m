function [] = Experiment0_Standard(ex_num, ex_num_variant, bool_preproc)
% 
%
% % Inputs
%
%
% ex_num : (String) Example number
%
% ex_num_variant : (String) 'a', 'b' or 'c' Changes the order of the three
%                   polynomials from the examples file.
%
% bool_preproc : (Boolean)
%       true : Include preprocessing of the subresutlant matrices
%       false : Exclude preprocessing of the subresutlant matrices
%
%
% % Examples
%
%
% >> Experiment0_Standard('1', 'a', true)
% >> Experiment0_Standard('1', 'b', true)
% >> Experiment0_Standard('1', 'c', true)

% Set upper and lower noise level
emin = 1e-10;
emax = 1e-12;


% Set preprocessing related variables
switch bool_preproc
    
    case true
        mean_method = "Geometric Mean Matlab Method";
        bool_alpha_theta = true;
        
    case false
        mean_method = "None";
        bool_alpha_theta = false;
        
end

% Low Rank Approximation Method
%   'None'
%   'Standard STLN'
%   'Standard SNTLN'
low_rank_approx_method = "None";

% APF Method
%   'None'
%   'Standard APF'
apf_method = "None";

% Sylvester matrix Variant
%   'T'
%   'DT'
%   'TQ'
%   'DTQ'
sylvester_matrix_variant = "DTQ";

% nEquations : Determines the structure of the three-polynomial subresulant
% matrix 
%       '2' : Subresultant matrices have a 2 by 3 partitioned structure
%       '3' : Subresutlant matrices have a 3 by 3 partitioned structure
nEquations = '2';

% Rank Revealing Metric
%   'Minimum Singular Values'
rank_revealing_metric = "Minimum Singular Values";


o_gcd_Univariate_3Polys(ex_num, ex_num_variant, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method,...
    sylvester_matrix_variant, nEquations, rank_revealing_metric)


end
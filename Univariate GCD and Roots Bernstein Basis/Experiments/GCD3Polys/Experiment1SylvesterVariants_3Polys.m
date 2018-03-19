function [] = Experiment1SylvesterVariants_3Polys(ex_num, ex_num_variant, bool_preproc)
% Experiment1SylvesterVariants_3Polys(ex_num, ex_num_variant, bool_preproc)
%
% % Inputs
%
%
% ex_num : (String) Example number
%
% ex_num_variant : (String) 'a', 'b' or 'c'
% 
% bool_preproc : (Boolean) 
%       true : Include preprocessing of the subresultant matrices
%       false : Exclude preprocessing of the subresultant matrices



% Set componentwise signal to noise upper and lower bounds
emin = 1e-8;
emax = 1e-8;

% Low Rank Approximation Method
%   'None'
%   'Standard STLN'
%   'Standard SNTLN'
low_rank_approx_method = 'None';


% Approximate Factorisation Method
%   'None'
%   'Standard APF'
apf_method = 'None';


% Rank Revealing Metric 
%   'Minimum Singular Values'
%   ''
rank_revealing_metric = 'Minimum Singular Values';


% nEquations : Used to determine the structure of the three-polynomial
% subresultant matrices
%       '2' : Subresultant matrices have a 2 by 3 partitioned structure.
%       '3' : Subresultant matrices have a 3 by 3 partitioned structure.
nEquations = '2';


% Set preprocessing related variables
switch bool_preproc
    
    case 1
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case 0
        mean_method = 'None';
        bool_alpha_theta = false;
    otherwise 
        error("Not  valid")
end



% % Variables tested in this file =========================================

% Initialise array of the Sylvester matrix variants
arrSylvesterMatrixVariants = {'T', 'DT', 'TQ', 'DTQ'};

% % =======================================================================





for i1 = 1 : 1 : length(arrSylvesterMatrixVariants)
    
    sylvester_matrix_variant = arrSylvesterMatrixVariants{i1};
    
    %try
    o_gcd_Univariate_3Polys(ex_num, ex_num_variant, emin, emax, mean_method, ...
        bool_alpha_theta, low_rank_approx_method, apf_method,...
        sylvester_matrix_variant, nEquations, rank_revealing_metric)
    %catch
    
    
    %end
    
    
end



end

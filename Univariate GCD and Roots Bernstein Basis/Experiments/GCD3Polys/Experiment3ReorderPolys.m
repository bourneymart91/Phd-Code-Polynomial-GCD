function [] = Experiment3ReorderPolys(ex_num, bool_preproc)
% This experiment file runs an example for three different orderings of the
% polynomials taken from the example file. For the first ordering 'a' 
% matches the polynomials (a(x),b(x),c(x)) to (f(x),g(x),h(x)). The second
% ordering 'b' matches (a(x),b(x),c(x)) to (g(x), f(x), h(x)). The third
% ordering 'c' matches (a(x),b(x),c(x)) to (h(x), g(x), f(x)). This is
% equivalent to considering the three variations of the (2 x 3) subresultant
% matrices.
% 
%
% % Inputs
%
%
% ex_num : (String) Example Number
%
% bool_preproc : (Boolean)
%
%
% % Outputs
%
% Results are printed to the .dat file called "Results_o_gcd_3Polys.dat"
% contained in the "Results" folder.
%
% % Examples 
%
% >> Experiment3ReorderPolys('1', true)


close all; 
clc;


% Set variables associated with preprocessing the subresultant matrices
%       true : Preprocess the subresultant matrices
%       false : Exclude preprocessing of the subresultant matrices

switch bool_preproc
    case true
        
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case false
        
        mean_method = 'None';
        bool_alpha_theta = false;
        
end


% Set upper and lower noise level
emin = 1e-12;
emax = 1e-10;


% Low Rank Approximation Method
%   'None'
%   'Standard STLN'
%   'Standard SNTLN'
low_rank_approx_method = 'None';


% APF Method
%   'None'
%   'Standard APF'
apf_method = 'None';


% nEquations
%   '2' : Sylvester subresultant matrices have a 2 by 3 structure
%   '3' : Sylvester subresultant matrices have a 3 by 3 structure
%
% Note - For this example, leave this set to '2' since we are experimenting
% with the optimal polynomial ordering in the (2 x 3) subresultant
% matrices. The (3 x 3) subresultant matrices are unaffected by reordered 
% polynomials.

nEquations = '2';



% Rank Revealing Metric
%   'Minimum Singular Values'
%   'R1 Row Diagonals'
rank_revealing_metric = 'Minimum Singular Values';



% Sylvester matrix variant
%   'T'
%   'DT'
%   'TQ'
%   'DTQ'
sylvester_matrix_variant = 'DTQ';

% Initialise array of example number variants
arr_ex_num_variant = {'a', 'b', 'c'};

for i = 1 : 1 : 3
    
    ex_num_variant = arr_ex_num_variant{i};
    
    o_gcd_Univariate_3Polys(ex_num, ex_num_variant, emin, emax, mean_method, ...
        bool_alpha_theta, low_rank_approx_method, apf_method,...
        sylvester_matrix_variant, nEquations, rank_revealing_metric)
end


end

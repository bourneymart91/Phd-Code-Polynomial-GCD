function [] = Experiment_As_Two_Poly_Problem(ex_num, bool_preproc)
% This experiment considers the GCD of three polynomials as two polynomial
% problems
% 
% % Inputs
%
%
% ex_num : (String) Example number
%
% bool_preproc : (Boolean) 
%       true : Include preprocessing of the subresultant matrices
%       false : Exclude preprocessing of the subresutlant matrices
%
% % Examples 
%
%
% >> Experiment_As_Two_Poly_Problem('11', false)




close all; clc;

% Set max and min noise level
emin = 1e-9;
emax = 1e-9;



% Set the preprocessing related variables
switch bool_preproc
    case true
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
    case false
        mean_method = 'None';
        bool_alpha_theta= false;
end

% Constants
apf_method = 'None';

% % Low Rank Approximation Method
%   'None'
%   'Standard STLN'
%   'Standard SNTLN'
low_rank_approx_method = 'None';

% Method used to determine the degree of the GCD
%   'R1 Row Norms',
%   'R1 Row Diagonals',
%   'Minimum Singular Values',
%   'Normalised Minimum Singular Values'
rank_revealing_metric = 'Minimum Singular Values';

% Set the sylvester matrix variant to be used. 
%   'T'
%   'DT'
%   'TQ'
%   'DTQ'
sylvester_matrix_variant = 'DTQ';

arr_ex_variant = {'a','b','c'};

for i = 1 : 1  : length(arr_ex_variant)
    
    ex_num_variant = arr_ex_variant{i};
    
    o_gcd_Univariate_2Polys(ex_num, ex_num_variant, emin, emax, mean_method, ...
        bool_alpha_theta, low_rank_approx_method, apf_method, ...
        sylvester_matrix_variant, rank_revealing_metric) ;
    
    
end






three_poly_problem = false;

if three_poly_problem == true
    
    % Three Polynomials - 2 x 3 subresutlant - f appears twice S(f,g,h)
    
    % nEquations refers to the choice of a (2x3) or a (3x3) partitioned
    % subresultant matrix.
    nEquations = '2';
    
    for i = 1 : 1 : length(arr_ex_variant)
        ex_variant = arr_ex_variant{i};
        ex_num = strcat(ex_num, ex_variant);
        o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, ...
            bool_alpha_theta, low_rank_approx_method, apf_method, ...
            sylvester_matrix_variant, nEquations, rank_revealing_metric)
    end
    
    
    
    % Three polynomials - 3 x 3 subresultant
    
    nEquations = '3';
    o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, ...
        bool_alpha_theta, low_rank_approx_method, apf_method, ...
        sylvester_matrix_variant, nEquations, rank_revealing_metric)
    
end




end
function [] = Experiment1SylvesterVariants_3Polys(ex_num)
% Experiment where alternate variatns of the subresultant matrices are
% considered.
%
% % Inputs
%
% ex_num : (String) Example number
%
%
% % Examples
%
%
% >> Experiment1SylvesterVariants_3Polys('14')



addpath(genpath(pwd));

%close all;
%clc;

% Good Examples
% ex_num = '14'; el = 1e-7; eu = 1e-7;
% ex_num = '14'; el = 1e-6; eu = 1e-6;

% Constants

% Set upper and lower limit
el = 1e-6;
eu = 1e-5;


bool_preprocessing = true;

switch bool_preprocessing
    case true
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case false
        
        mean_method = 'None';
        bool_alpha_theta = false;
end

% Low Rank Approximation Method
% 'None'
% 'Standard STLN'
low_rank_approx_method = 'None';

% Approximate Polynomial Factorisation Method
% 'None'
apf_method = 'None';

% Rank Revealing Metric
% 'Minimum Singular Values'
rank_revealing_metric = 'Minimum Singular Values';


%
arrSylvesterFormat = {...
    'T', ...
    'DT', ...
    'TQ', ...
    'DTQ'};

% nEquations determines the 
nEquations = '3';

for i = 1 : 1 : length(arrSylvesterFormat)
    
    % Get variant of subresultant matrix
    sylvester_format = arrSylvesterFormat{i};
    
    % Compute the GCD
    o_gcd_Bivariate_3Polys(ex_num, el, eu, mean_method, bool_alpha_theta, ...
        low_rank_approx_method, apf_method, sylvester_format, ...
        rank_revealing_metric, nEquations)
    
end

end
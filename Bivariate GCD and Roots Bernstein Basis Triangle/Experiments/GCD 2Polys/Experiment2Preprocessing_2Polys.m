function [] = Experiment2Preprocessing_2Polys(ex_num)
% Experiment computing the GCD with and without preprocessing the
% polynomials f(x,y) and g(x,y)
%
% % Inputs
%
%
% ex_num : (String) Example Number
%
%
% % Examples
%
%
% >> Experiment2Preprocessing_2Polys('19')



close all;
clc;

% Good Examples



% Constants

% Set upper and lower noise level
el = 1e-6;
eu = 1e-5;

% Low rank approximation method
% 'None'
% 'Standard STLN'
low_rank_approx_method = 'None';

%
apf_method = 'None';

% Rank Revealing Metric
% 'Minimum Singular Values'
% ''
rank_revealing_metric = 'Minimum Singular Values';

% Sylvester Format
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_format = 'DTQ';


% Without Preprocessing

mean_method = 'None';
bool_alpha_theta = false;

o_gcd_Bivariate_2Polys(ex_num, el, eu, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric)



% With Preprocessing

mean_method = 'Geometric Mean Matlab Method';
bool_alpha_theta = true;

o_gcd_Bivariate_2Polys(ex_num, el, eu, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric)


end

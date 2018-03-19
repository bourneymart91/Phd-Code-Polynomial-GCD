function [] = Experiment2Preprocessing_3Polys(ex_num)
%
%
% Experiment2Preprocessing_3Polys('7')

addpath(genpath(pwd));

close all; 
clc;

%
% 1 :
% 2 :
% 3 : Univariate GCD
% 4 : Univariate GCD
% 5 : Too Small
%
%
%

% Good Examples
% ex_num = '7'; el = 1e-4; eu = 1e-6;

% ex_num = '13'; el = 1e-6; eu = 1e-8;

% Constants
el = 1e-5;
eu = 1e-5;

low_rank_approx_method = 'None';
apf_method = 'None';
rank_revealing_metric = 'Minimum Singular Values';

sylvester_format = 'DTQ';
nEquations = '2';


mean_method = 'None';
bool_preproc = false;

o_gcd_Bivariate_3Polys(ex_num, el, eu, mean_method, bool_preproc, ...
    low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric, nEquations)

mean_method = 'Geometric Mean Matlab Method';
bool_preproc = true;

o_gcd_Bivariate_3Polys(ex_num, el, eu, mean_method, bool_preproc, ...
    low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric, nEquations)


end
function [] = Experiment3LowRankApproximation(ex_num)

addpath(genpath(pwd));

close all; 
clc;

% Constants
el = 1e-8;
eu = 1e-8;

apf_method = 'None';
rank_revealing_metric = 'Minimum Singular Values';

sylvester_format = 'DTQ';

mean_method = 'Geometric Mean Matlab Method';
bool_preproc = true;


low_rank_approx_method = 'None';
o_gcd_Bivariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric)

low_rank_approx_method = 'Standard STLN';
o_gcd_Bivariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric)

%low_rank_approx_method = 'Standard SNTLN';
%o_gcd_Bivariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric)



end
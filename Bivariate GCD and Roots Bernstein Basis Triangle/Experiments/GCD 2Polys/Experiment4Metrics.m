addpath(genpath(pwd));

close all; 
clc;

% Constants
ex_num = '12';
el = 1e-10;
eu = 1e-10;

low_rank_approx_method = 'None';
apf_method = 'None';

sylvester_format = 'DTQ';



%mean_method = 'None';
%bool_preproc = false;
mean_method = 'Geometric Mean Matlab Method';
bool_preproc = true;




rank_revealing_metric = 'Minimum Singular Values';
o_gcd_Bivariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric)


rank_revealing_metric = 'Max Min R1 Row Diagonals';
o_gcd_Bivariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric)

rank_revealing_metric = 'Residuals';
o_gcd_Bivariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric)


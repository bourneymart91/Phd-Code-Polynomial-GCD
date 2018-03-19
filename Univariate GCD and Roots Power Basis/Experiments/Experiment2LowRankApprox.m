
close all; clc;


ex_num = '6';
emin = 0;
emax = 0;
mean_method = 'Geometric Mean Matlab Method';
bool_alpha_theta = true;
%mean_method = 'None';
%bool_alpha_theta = false;
low_rank_approx_method = 'Standard STLN';
rank_revealing_metric = 'Minimum Singular Values';


apf_method = 'None';
o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, rank_revealing_metric)

apf_method = 'Last Non-zero Row';
o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, rank_revealing_metric)
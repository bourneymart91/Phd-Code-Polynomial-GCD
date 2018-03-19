function [] = Experiment2_3Polys(ex_num, bool_preproc)
close all; clc;

el = 1e-10;
eu = 1e-8;

if bool_preproc == true
    
    mean_method = 'Geometric Mean Matlab Method';
    bool_alpha_theta = true;
else
    mean_method = 'None';
    bool_alpha_theta = false;
end


subresultant_build_method = 'T';
rank_revealing_metric = 'Minimum Singular Values';
low_rank_approx_method = 'None';
apf_method = 'None';



% Compute the GCD using only the (2 x 3) partitioned subresultant matrices
nEquations = '2';

ex_num_variation = strcat(ex_num,'a');
o_gcd_Univariate_3Polys(ex_num_variation, el, eu, mean_method, bool_alpha_theta, subresultant_build_method , low_rank_approx_method , apf_method, rank_revealing_metric, nEquations);

ex_num_variation = strcat(ex_num,'b');
o_gcd_Univariate_3Polys(ex_num_variation, el, eu, mean_method, bool_alpha_theta, subresultant_build_method , low_rank_approx_method , apf_method, rank_revealing_metric, nEquations);

ex_num_variation = strcat(ex_num,'c');
o_gcd_Univariate_3Polys(ex_num_variation, el, eu, mean_method, bool_alpha_theta, subresultant_build_method , low_rank_approx_method , apf_method, rank_revealing_metric, nEquations);

% Compute the GCD using only the (3 x 3) partitioned subresultant matrices
nEquations = '3';

ex_num_variation = strcat(ex_num,'a');
o_gcd_Univariate_3Polys(ex_num_variation, el, eu, mean_method, bool_alpha_theta, subresultant_build_method , low_rank_approx_method , apf_method, rank_revealing_metric, nEquations);



end
function [] = Experiment_As_Two_Poly_Problem(ex_number, bool_preproc)

% Variable : Preprocessing

% Good Examples

% Example 11 - 1e-10 , 1e-12



close all; clc;

emin = 1e-9;
emax = 1e-9;

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
low_rank_approx_method = 'None';
rank_revealing_metric = 'Minimum Singular Values';
sylvester_matrix_variant = 'DTQ';

arr_ex_variant = {'a','b','c'};

for i = 1 : 1  : length(arr_ex_variant)
    
    ex_variant = arr_ex_variant{i}; 
    ex_num = strcat(ex_number, ex_variant);
    o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_matrix_variant, rank_revealing_metric) ;


end

nEquations = '2';

% Three Polynomials - 2 x 3 subresutlant - f appears twice S(f,g,h)
for i = 1 : 1 : length(arr_ex_variant)
    ex_variant = arr_ex_variant{i}; 
    ex_num = strcat(ex_number, ex_variant);
    o_gcd_Bivariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_matrix_variant, rank_revealing_metric, nEquations)
end



% Three polynomials - 3 x 3 subresultant 

nEquations = '3';
o_gcd_Bivariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_matrix_variant, rank_revealing_metric, nEquations)






end
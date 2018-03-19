function [] = Experiment3_PolyOrdering(ex_number, bool_preproc)
%
%
% Experiment3PolyOrdering

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



switch bool_preproc
    
    case true
        mean_method = 'Geometric Mean Matlab Method';
        bool_preproc = true;
        
    case false
        mean_method = 'None';
        bool_preproc = false;
end

% Constants
el = 1e-8;
eu = 1e-6;
low_rank_approx_method = 'None';
apf_method = 'None';
rank_revealing_metric = 'Minimum Singular Values';
sylvester_format = 'DTQ';


arrVariants = {'a','b','c'};
nVariants = length(arrVariants);


nEquations = '2';
for i = 1 : 1 : nVariants
   
    variant = arrVariants{i};
    ex_num = strcat(ex_number, variant);
    o_gcd_Bivariate_3Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric, nEquations)

end


nEquations = '3';
o_gcd_Bivariate_3Polys(ex_num, el, eu, mean_method, bool_preproc, low_rank_approx_method, apf_method, sylvester_format, rank_revealing_metric, nEquations)





end
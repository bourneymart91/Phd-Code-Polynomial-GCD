function [] = Experiment7_DegreeElevation(ex_num, bool_preproc, p, q)
% This experiment looks at the effect of degree elevation in the
% computation of the degree of the GCD of two univariate polynomials in
% Bernstein form.
%
%
% % Inputs
%
%
% ex_num : (String) Example number
%
% bool_preproc : (Boolean) 
%
% p : (Int) Number of degree elevations of f(x)
%
% q : (Int) Number of degree elevations of g(x)
%
% 
%
%
% % Examples
%
% >> Experiment7_DegreeElevation('1', true, 5, 2)




% Constants --------------------------------------------------------------


% Set min and max noise level
emin = 1e-12;
emax = 1e-12;

% Set preprocessing related variables
switch bool_preproc
    case true      
        
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case false
        
        mean_method = 'None';
        bool_alpha_theta = false;
        
end

% Set other variables
low_rank_approx_method = 'None';
apf_method = 'None';

% Set the sylvester matrix variant to be used. 
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_matrix_variant = 'DTQ';


rank_revealing_metric = 'Minimum Singular Values';

% Compute the GCD 
o_gcd_Univariate_2Polys_DegreeElevationTest(ex_num, emin, emax, ...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_variant, rank_revealing_metric, p, q)

end
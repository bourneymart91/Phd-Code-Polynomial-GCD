function [] = o_roots_Bivariate(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, degree_method, ...
    rank_revealing_metric, deconvolution_method_hxy, deconvolution_method_wxy)
% o_roots_bivar(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
%
% Given an example number and set of parameters, obtain the roots of the
% example polynomial, where the polynomial is in the power basis
%
% Inputs
%
% ex_num : (String) Example Number
%
% emin : (Float) Lower noise level
%
% emax : (Float) Upper noise level
%
% mean_method : 
%       'Geometric Mean Matlab Method'
%       'None'
%
% bool_alpha_theta (Bool)
%       true - Include Preprocessing
%       false - Exclude Preprocessing
%
% low_rank_approx_method 
%       'Standard SNTLN' : Include Nonlinear SNTLN
%       'Standard STLN'  : Standard Linear STLN
%       'None'           : Exclude SNTLN
%
% apf_method
%       'Standard APF'
%       'None'
%
% degree_method
%       'Total'
%       'Relative'
%       'Both'
%
% rank_revealing_metric
%
% % Examples
% >> o_roots_Bivariate('1', 1e-10, 1e-12, 'None', false, 'None', 'None', 'Both', 'Minimum Singular Values', 'Batch', 'Separate')
% >> o_roots_Bivariate('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true , 'Standard STLN', 'Standard APF', 'Both', 'Minimum Singular Values', 'Batch', 'Separate')

%
restoredefaultpath
addpath(genpath(pwd));



if emax < emin
   temp_min = emax;
   emax = emin;
   emin = temp_min;
   
end

SetGlobalVariables_Roots( ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, degree_method, ...
    rank_revealing_metric, deconvolution_method_hxy, deconvolution_method_wxy);

% %
% %
% Get Example

% Given the example number get the polynomial coefficients of f(x,y) and
% the total degree 'm' of f(x,y).
[fxy_exact, m] = Examples_Roots_Bivariate(ex_num);

% Add noise to the coefficients of polynomial f(x,y)
[fxy,~] = AddVariableNoiseToPoly(fxy_exact, emin, emax);


% % 
% %
% Get roots by my method

root_finding_method = '3 Poly GCD';

switch root_finding_method 
    case '2 Poly GCD'
        
        o_roots_mymethod_Bivariate(fxy, m);        
        
    case '3 Poly GCD'

        o_roots_mymethod_newmethod_Bivariate(fxy, m);
        
    otherwise
        
        error([mfilename ' : Error \n'])
end

end
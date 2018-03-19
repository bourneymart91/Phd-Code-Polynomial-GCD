function [] = SetGlobalVariables_GCD_3Polys( ex_num, ex_num_variant, emin, emax, ...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method,...
    sylvester_matrix_variant, nEquations, rank_revealing_metric)
% Set the global variables
%
% Inputs.
%
% ex_num : (String) Example number
%
% ex_num_variant : (String) 'a', 'b' or 'c' determines ordering of the
% polynomials f(x), g(x) and h(x)
%
% emin : (Float) Mininimum componentwise signal to noise
%
% emax : (Float) Maximum componentwise signal to noise
%
% mean_method : (String)
%   * 'None'
%   * 'Geometric Mean Matlab Method'
% 
% bool_alpha_theta : (Boolean)
%   * true
%   * false
%
% low_rank_approx_method : (String)
%   * None :   No perturbations added
%   * Standard STLN :   Use linear form
%   * Root Specific SNTLN :   Non-Linear low form, with g = f' constraint
%   * Standard SNTLN :   Non-Linear low rank form, without constraint
%
% apf_method : (String)
%
% Sylvester_matrix_variant : (String)
%   * T : 
%   * DT :
%   * DTQ :
%   * TQ : 
%   * DTQ Denominator Removed :
%   * DTQ Rearranged :
%
% nEquaitons : (String)
% number of equations used to define the sylvester matrix
%   * '2' : Sylvester matrix is 2 x 3 structure
%   * '3' : Sylvester matrix is 3 x 3 structure
%
% rank_revealing_metric : (String)
%   * Singular Values :   
%   * Max/Min Singular Values :
%   * R1 Row Norms : 
%   * R1 Row Diagonals :
%   * Residuals :

if (nargin ~= 11)
   error('Custom Alert : Not enough input arguments') 
end

global SETTINGS

SETTINGS.EX_NUM = ex_num;
SETTINGS.EX_NUM_VARIANT = ex_num_variant;
SETTINGS.EMIN = emin;
SETTINGS.EMAX = emax;
SETTINGS.SEED = 1024;


% PLOT_GRAPHS
SETTINGS.PLOT_GRAPHS = true;
SETTINGS.PLOT_GRAPHS_GCD_DEGREE = true;
SETTINGS.PLOT_GRAPHS_PREPROCESSING = false;
SETTINGS.PLOT_GRAPHS_LOW_RANK_APPROXIMATION = false;

% BOOL_LOG
%   * true : Use logs
%   * false : Dont use logs
%
SETTINGS.BOOL_LOG = false;

%--------------------------------------------------------------------------
%               
% Preprocessing related values
%
%
SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;
SETTINGS.SCALING_METHOD = 'lambda_mu_rho';
SETTINGS.MEAN_METHOD = mean_method;

% -------------------------------------------------------------------------
% 
% SETTINGS IN COMPUTING GCD DEGREE
%

% Set the metric for measuring the degree of the GCD.
%
%   * Singular Values :   
%   * Max/Min Singular Values :
%   * R1 Row Norms : 
%   * R1 Row Diagonals :
%   * Residuals :
%
SETTINGS.RANK_REVEALING_METRIC = rank_revealing_metric;


% -------------------------------------------------------------------------
%
% SETTINGS FOR COMPUTING GCD COEFFICIENTS
%
%

% When computing the GCD do we wish to use both u(x,y) and v(x,y) or only
% u(x,y). Note: It is necessary to use both for standard GCD computations.
%
% 'ux and vx' :
% 'ux' :
%
SETTINGS.GCD_COEFFICIENT_METHOD = 'ux and vx';


%--------------------------------------------------------------------------

% Structuring the Sylvester Matrix
% SYLVESTER_MATRIX_VARIANT
%   * T : 
%   * DT :
%   * DTQ :
%   * TQ : 
%   * DTQ Denominator Removed :
%   * DTQ Rearranged :
SETTINGS.SYLVESTER_MATRIX_VARIANT = sylvester_matrix_variant;


% SYLVESTER_EQUATIONS 
%   * '2' : Sylvester Matrix defined by 2 equations
%   * '3' : Sylvester Matrix defined by 3 equations
SETTINGS.SYLVESTER_EQUATIONS = nEquations;

%-------------------------------------------------------------------------
%
% LOW_RANK_APPROXIMATION_METHOD:
%
%   * None : No perturbations added
%   * Standard STLN : Use linear form
%   * Root Specific SNTLN : Non-Linear low form, with g = f' constraint
%   * Standard SNTLN : Non-Linear low rank form, without constraint
%
SETTINGS.LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;

SETTINGS.MAX_ERROR_SNTLN = 1e-13;
SETTINGS.MAX_ITERATIONS_SNTLN = 20;

%--------------------------------------------------------------------------
%
% APF RELATED SETTINGS.
%
%

% Structuring the matrix [C(f) | C(g)]
%   * 'Standard APF NonLinear'
%   * 'Standard APF Linear'
%   * 'None'
SETTINGS.APF_METHOD = apf_method;

% APF_BUILD_METHOD
% * Standard
% * Rearranged
SETTINGS.APF_BUILD_METHOD = 'Standard';

% Regarding the computation of the low rank approximation of the 
% C = [C(u) ; C(v)] matrix
SETTINGS.MAX_ERROR_APF = 1e-14;
SETTINGS.MAX_ITERATIONS_APF = 50;





% -------------------------------------------------------------------------
% Print the parameters to console
PrintParameters_GCD_3Polys()

end





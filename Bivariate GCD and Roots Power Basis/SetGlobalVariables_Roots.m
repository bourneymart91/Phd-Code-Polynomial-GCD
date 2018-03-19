function [] =  SetGlobalVariables_Roots(problem_type, ex_num, emin, emax, ...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, ...
    degree_method, rank_revealing_metric, deconvolution_method_hxy, ...
    deconvolution_method_wxy)
% SetGlobalVariables()
%
%
% % Inputs
%
% problem_type : (String) Either GCD or Roots
%
% ex_num : (String)
%
% emin : (Float)
%
% emax : (Float) 
%
% mean_method : (String)
%   * Geometric Mean Matlab Method
%   * None
%
% bool_alpha_theta : (Boolean)
%   * true
%   * false
%
%
% low_rank_approx_method : (String)
%   * None
%   * Standard STLN
%   * Standard SNTLN
%
% apf_method
%   * None
%   * Linear APF
%   * Nonlinear APF
%
% degree_method
%   * Total
%   * Relative
%   * Both
%
% rank_revealing_metric
%   * R1 Row Norms
%   * R1 Row Diagonals
%   * Singular Values
%   * Residuals
%
%
% deconvolution_method_hxy
%  * Separate
%  * Batch
%
%
% deconvolution_method_wxy
%  * Separate
%  * Batch

global SETTINGS

SETTINGS.PLOT_GRAPHS = true;

% Problem Type : 'GCD' or 'Roots'
SETTINGS.PROBLEM_TYPE = problem_type;

% SETTINGS.VECTORISATION METHOD
%   Version 1 - Read diagonally
%   Version 2 - Read down rows then across columns
SETTINGS.VECTORISATION_METHOD = 'Version 2';

%-------------------------------------------------------------------------
%
% Example Settings
%
%

% Example Number : String
SETTINGS.EX_NUM = ex_num;

% Noise upper and lower limit
SETTINGS.EMIN = emin;
SETTINGS.EMAX = emax;

SETTINGS.SEED = 1024;

%--------------------------------------------------------------------------


% Degree_method
%       Total
%       Relative
%       Both

SETTINGS.DEGREE_METHOD = degree_method;

%--------------------------------------------------------------------------
%
%       SETTINGS : PREPROCESSING
%
%


SETTINGS.MEAN_METHOD = mean_method;

SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;


%-------------------------------------------------------------------------
%
%       SETTINGS : DEGREE COMPUTATION
%
%
SETTINGS.THRESHOLD = 2;
SETTINGS.THRESHOLD_RANK = 2;


% Metric used to compute the degree of the GCD
% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals

SETTINGS.RANK_REVEALING_METRIC = rank_revealing_metric;

%--------------------------------------------------------------------------
%
% % SETTINGS : LOW RANK APPROXIMATION
%
%

% Low Rank Approximation Method
%       'Standard STLN'
%       'Standard SNTLN'
%       'None'
SETTINGS.LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;

SETTINGS.MAX_ITERATIONS_SNTLN = 50;

SETTINGS.MAX_ERROR_SNTLN = 1e-10;

%--------------------------------------------------------------------------
%
% SETTINGS : APPROXIMATE FACTORISATION SETTINGS
%
%
SETTINGS.APF_METHOD = apf_method;

%--------------------------------------------------------------------------
%
% SETTINGS : DECONVOLUTION
% Batch - Perform batch of deconvolutions together
% Separate - Each deconvolution is separate
%
SETTINGS.DECONVOLUTION_METHOD_HXY = deconvolution_method_hxy;

%
% 'From Deconvolutions'
% 'From uxy'
%
SETTINGS.HXY_METHOD = 'From uxy';
SETTINGS.DECONVOLUTION_METHOD_WXY = deconvolution_method_wxy;

end


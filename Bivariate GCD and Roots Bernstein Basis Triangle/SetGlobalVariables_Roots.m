function [] = SetGlobalVariables_Roots( ex_num, emin, emax,...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_variant, rank_revealing_metric, ...
    deconvolution_method_hx, deconvolution_method_wx, nEquations)
%
% % Inputs
%
% ex_num : (String) Example Number
%
% emin : (Float)
%
% emax : (Float)
%
% mean_method : (String)
%   * None
%   * Geometric Mean Matlab Method
%
% bool_alpha_theta : (Boolean)
%   * true
%   * false
%
% low_rank_approx_method : (String)
%   * 'Standard STLN'
%   * 'None'
%
% apf_method : (String)
%
% sylvester_matrix_variant : (String)
%   * R1 Row Norms
%   * R1 Row Diagonals
%   * Singular Values
%   * Residuals
%
%
% deconvolution_method_hx : (String)
%   * Separate
%   * Batch
%   * Batch With STLN
%   * Batch Constrained
%   * Batch Constrained With STLN
%
%
% deconvolution_method_wx : (String)
%   * Separate
%   * Batch
%   * Batch With STLN
 


global SETTINGS

%----------------------

% Plot Graphs
%   * true
%   * false
%
SETTINGS.PLOT_GRAPHS = true;



SETTINGS.SEED = 1024;

% SYLVESTER BUILD METHOD
%   T
%   DTQ
%   DT
%   TQ
%   DTQ Denominator Removed

SETTINGS.SYLVESTER_MATRIX_VARIANT = sylvester_matrix_variant;
SETTINGS.SYLVESTER_MATRIX_3POLY_N_EQUATIONS = nEquations;

% Set Noise
SETTINGS.EMIN = emin;
SETTINGS.EMAX = emax;

% Set example number
SETTINGS.EX_NUM = ex_num;

%--------------------------------------------------------------------------
%
%       SETTINGS - PREPROCESSING
%
%


% MEAN_METHOD
% 'Geometric Mean Matlab Method'
% 'Geometric Mean My Method'
% 'None'
%
SETTINGS.MEAN_METHOD = mean_method;

% BOOL_ALPHA_THETA
%   * true
%   * false
%
SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;

% lambda_mu_rho
% lambda_mu
% lambda_rho
% mu_rho
SETTINGS.SCALING_METHOD = 'lambda_rho';

% -------------------------------------------------------------------------
%
%   SETTINGS - DEGREE COMPUTATION
%
%
SETTINGS.THRESHOLD = 2;
SETTINGS.THRESHOLD_RANK = 1e-3;

% RANK_REVEALING_METRIC
%   * R1 Row Norms
%   * R1 Row Diagonals
%   * Singular Values
%   * Residuals
SETTINGS.RANK_REVEALING_METRIC = rank_revealing_metric;

% ------------------------------------------------------------------------
%
%   SETTINGS - LOW RANK APPROXIMATION
%
% LOW_RANK_APPROX_METHOD
% 'Standard STLN'
% 'None'
%
SETTINGS.LOW_RANK_APPROX_METHOD = low_rank_approx_method;

SETTINGS.STLN_MAX_ERROR = 1e-10;
SETTINGS.STLN_MAX_ITERATIONS = 50;

%--------------------------------------------------------------------------
%
%   SETTINGS - APPROXIMATE FACTORISATION
%
%
SETTINGS.APF_METHOD = apf_method;


%-------------------------------------------------------------------------
%
%   SETTINGS - DECONVOLUTION
%
%
%

% Deconvolution Settings for root finding method

% Method used to deconvolve polynomials f_{i}(x,y) to compute h(x,y)
%
% DECONVOLUTION_METHOD_HXY
%   * 'Separate'
%   * 'Batch'
%   * 'Batch With STLN'
%   * 'Batch Constrained'
%   * 'Batch Constrained With STLN'
%
SETTINGS.DECONVOLUTION_METHOD_HXY = deconvolution_method_hx;


% Method used to deconvolve polynomials h_{i}(x,y) to compute w_{i}(x,y)
%
%   * 'Separate'
%   * 'Batch'
%   * 'Batch With STLN'
%
SETTINGS.DECONVOLUTION_METHOD_WXY = deconvolution_method_wx;
SETTINGS.PREPROC_DECONVOLUTIONS = true;

SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-10;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 20;
end

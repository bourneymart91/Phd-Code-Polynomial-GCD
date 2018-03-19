function [] = SetGlobalVariables_GCD_2Polys( ex_num, emin, emax,...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_variant, rank_revealing_metric)
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
%
% bool_alpha_theta : (Boolean)
%   * True
%   * False
%
% low_rank_approx_method : (String)
%   * None
%   * Standard STLN
%   * Standard SNTLN
%
%
% apf_method : (String)
%
% sylvester_matrix_variant : (String)
%   * T
%   * DT
%   * TQ
%   * DTQ
%   * DTQ Denominator Removed
%
% rank_revealing_metric : (String)
%   
%

global SETTINGS

%----------------------

% Plot Graphs
%   * true
%   * false
%
SETTINGS.PLOT_GRAPHS = true;
SETTINGS.PLOT_GRAPHS_PREPROCESSING = false;
SETTINGS.PLOT_GRAPHS_DECONVOLUTION = false;
SETTINGS.PLOT_GRAPHS_RANK = true;
SETTINGS.PLOT_GRAPHS_LRA = false;

SETTINGS.SEED = 1024;

% SYLVESTER BUILD METHOD
%   T
%   DTQ
%   DT
%   TQ
%   DTQ Denominator Removed

SETTINGS.SYLVESTER_MATRIX_VARIANT = sylvester_matrix_variant;

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


%
% 'Geometric Mean Matlab Method'
% 'Geometric Mean My Method'
% 'None'
%
SETTINGS.MEAN_METHOD = mean_method;

%
%   * true
%   * false
%
SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;



% -------------------------------------------------------------------------
%
%   SETTINGS - DEGREE COMPUTATION
%
%
SETTINGS.THRESHOLD = 2;
SETTINGS.THRESHOLD_RANK = 1e-3;

% Metric used to compute the degree of the GCD
% R1 Row Norms
% R1 Row Diagonals
% Minimum Singular Values
% Residuals
SETTINGS.RANK_REVEALING_METRIC = rank_revealing_metric;

% ------------------------------------------------------------------------
%
%   SETTINGS - LOW RANK APPROXIMATION
%
%


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



end

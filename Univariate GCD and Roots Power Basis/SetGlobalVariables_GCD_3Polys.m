function [] = SetGlobalVariables_GCDFinding(ex_num, emin, emax, ...
    mean_method, bool_alpha_theta, sylvester_build_method, ...
    low_rank_approx_method, apf_method, rank_revealing_metric, ...
    nEquations)
%
%
% % Inputs
%
% ex_num : (String)
%
% emin : (Float)
%
% emax : (Float)
%
% mean_method : (String)
%
% bool_alpha_theta : (Boolean)
%
% sylvester_build_method : (String)
%
% low_rank_approx_method : (String)
%
% apf_method : (String)
%
% rank_revealing_metric : (String)
%
% nEquations
% 



global SETTINGS


% Example Number
SETTINGS.EX_NUM = ex_num;

% Seed for noise generation
SETTINGS.SEED = 1024;

% Minimum/Maximum noise level
SETTINGS.EMIN = emin;
SETTINGS.EMAX = emax;

% Include/Exclude plotting of graphs
SETTINGS.PLOT_GRAPHS = true;


SETTINGS.GCD_COEFFICIENT_METHOD = 'ux and vx';

%--------------------------------------------------------------------------
%
%           SETTINGS : PREPROCESSING
%
%

% Method of computing mean
SETTINGS.MEAN_METHOD = mean_method;

% Include/Exclude preprocessing
SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;

SETTINGS.SYLVESTER_BUILD_METHOD = sylvester_build_method;
SETTINGS.N_EQUATIONS_SYLVESTER_MATRIX = nEquations;
%--------------------------------------------------------------------------
%
%       SETTINGS : DEGREE COMPUTATION
%
%
%
%


% Threshold defined for use in computing the degree of the GCD
SETTINGS.THRESHOLD = 2;
SETTINGS.THRESHOLD_RANK = -5;

% (Roots Only)
% Make use of precalculated limits on the upper and lower bounds of the
% degree t of the GCD(f,f').
SETTINGS.BOOL_LIMITS = 'y';

% Metric used to compute the degree of the GCD
%   * Minimum Singular Values
%   * R1 Row Norms
%   * R1 Row Diagonals
%   * Residuals

SETTINGS.METRIC = rank_revealing_metric;

%--------------------------------------------------------------------------
%
%       SETTINGS : Low rank approximation method
%
%
%

%
% STLN
% SNTLN
SETTINGS.LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;

% Settings for SNTLN/STLN
SETTINGS.MAX_ERROR_SNTLN = 1e-12;
SETTINGS.MAX_ITE_SNTLN = 50;

%--------------------------------------------------------------------------
%
%       SETTINGS : APF Method
%
%
%

SETTINGS.APF_METHOD = apf_method;


end
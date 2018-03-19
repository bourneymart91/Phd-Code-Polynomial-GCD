function [] = Experiment_DeconvolutionMethod(ex_num)
% Experiment on the differente deconvolution method, with and without
% preprocessing
%
% % Inputs
%
% ex_num : (String) Example Number


close all; 
clc; 


global SETTINGS
SETTINGS.PLOT_GRAPHS_DECONVOLUTION = false;

% % Good Examples
%
% 1 : Small example, Great results for batch methods (1e-8) 
%
% 2 : Good example, shows good effect of preproc (1e-6)
%
% 3 : 
%
% 4 : 
%
% % Bad Examples
% 
% 5
%
%
%

% Set noise level
emin = 1e-6;

% Perform deconvolutions without preprocessing
bool_preproc = false;            
o_deconvolution(ex_num, emin, bool_preproc)

% Perform deconvolutions with preprocessing
bool_preproc = true;            
o_deconvolution(ex_num, emin, bool_preproc)



end
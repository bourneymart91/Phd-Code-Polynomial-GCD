function [] = Experiment1Preprocessing(ex_num)
% Perform a deconvolution problem of the type found in the root finding 
% algorithm. Given an example number, perform the deconvolution problem 
% with and without preprocessing. 
%
% % Inputs
%
% ex_num : (String) Example number.
%
%
% % Examples
%
% Experiment1Preprocessing('1')



close all;
clc;

% % Set noise level
noise = 1e-4;



% Test 1 - Perform deconvolution without preprocessing
bool_preproc = false;
o_Deconvolution(ex_num, noise, bool_preproc)

% Test 2 - Perform deconvolution with preprocessing
bool_preproc = true;
o_Deconvolution(ex_num, noise, bool_preproc)

end
function [] = Experiment2VariableNoise(ex_num, bool_preproc)
% Perform a series of deconvolutions as found in the factorisation
% algorithm, for a variety of noise levels.
%
% % Inputs
%
% ex_num : (String) Example number
%
% bool_preproc : (boolean) Set whether to include or exclude preprocessing 
%                   of the polynomials {f_{i}(x)}
%
%
% % Outputs
%
% Outputs are saved to a .dat file


close all;
clc;


% Create an array of noise levels, where noise is added to the coefficients
% of the set of polynomials f_{i}(x)
arrNoise = {1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4};



% For each noise level, perform the deconvolution
for i = 1 : 1 : length(arrNoise)
    
    noise = arrNoise{i};
    
    o_Deconvolution(ex_num, noise, bool_preproc)
    
    
end

end
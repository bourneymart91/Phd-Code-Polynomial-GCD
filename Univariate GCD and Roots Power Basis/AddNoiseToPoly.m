function [f_noisy, noise_vec] = AddNoiseToPoly(fx, emin)
% Given a polynomial f(x), add noise to its coefficients at the level
% specified by mue.
%
% % Inputs;
%
% fx : (Vector) Column vector of coefficients of polynomial f(x).
%
% emin : Lower noise level.
%
% % Outputs
%
% f_noisy : Column vector of noisy coefficients of polynomial f(x)
%
% noise_vec : Noise added to coefficients of f(x)
%
%
% % Example
% 
% Add noise to a polynomial at a level of 1e-10.
% >> Noise([1.1 ; 2.2 ; 3.3], 1e-10)

% Initialise the global variables
global SETTINGS

% Set the random number generator.
rng(SETTINGS.SEED)

% Get the degree of the polynomial f(x)
m = GetDegree(fx);

% Generate a vector of random numbers between -1 and 1
a = -1;
b = 1;
r = a + (b-a).*rand(m+1,1);

% Get the vector of noise
noise_vec = fx.*(r.*emin);

% Add noise to exact coefficients
f_noisy = fx + noise_vec;


end
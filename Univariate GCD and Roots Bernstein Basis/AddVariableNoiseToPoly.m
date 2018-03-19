function [f_noisy,noisevector] = AddVariableNoiseToPoly(fx, el, eu)
% This function adds noise in the componentwise sense to the coefficients
% of the Bernstein basis polynomial.
% The upper threshold \epsilon is a random variable between two values
% \epsilon_{max} and \epsilon_{min}.
%
%
% Inputs:
%
%
% fx : (Vector) exact coefficients of polynomial f, in the Bernstein basis.
%
% el : (Float) Lower threshold of the epsilon value (Noise/Signal)
%
% eu : (Float) Upper threshold of the epsilon value (Noise/Signal)
%
%
% Outputs:
%
%
% f_noisy : (Vector) The noisy coefficients of perturbed polynomial f.
%
% noisevector : (Vector) vector of the noise added to f_exact.


global SETTINGS

% Get degree of input polynomial f.
m = GetDegree(fx);

% Set Seed for random number generator.
rng(SETTINGS.SEED)

% Generate random variables r_{i}
r = (2*rand(m+1,1))-ones(m+1,1);


% Generate random variables v_{i}
v = rand(m+1, 1);

% Get vector 'eps' - the vector which stores the upper Noise/Signal
% threshold \epsilon_{i} for each cofficient a_{i}
eps = el + v.*(eu-el);

% Calculate the noise vector = a_{i}r_{i}.*\epsilon_{i}
noisevector = fx.* r .*eps;

% Calculate the perturbed coefficients = a_{i} + a_{i}r_{i}\epsilon_{i}
f_noisy = fx + noisevector;




end
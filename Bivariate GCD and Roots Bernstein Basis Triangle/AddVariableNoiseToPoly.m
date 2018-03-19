function [f_noisy, noise_matrix] = AddVariableNoiseToPoly(fxy, el, eu)
% This function adds noise in the componentwise sense to the coefficients
% of the Bernstein basis polynomial.
% The upper threshold \epsilon is a random variable between two values
% \epsilon_{max} and \epsilon_{min}.
%
% Inputs:
%
%
% fxy : (Matrix) Coefficients of the polynomial exact polynomial f(x,y).
%
% el : (Float) Lower threshold of the epsilon value (Noise/Signal)
%
% eu : (Float) Upper threshold of the epsilon value (Noise/Signal)
%
%
% Outputs:
%
%
% f_noisy : (Matrix) Noisy coefficients of perturbed polynomial f.
%
% noise_matrix : (Matrix)  the noise added to f_exact.


global SETTINGS

% Get degree of input polynomial f(x,y)
[m,~] = GetDegree_Bivariate(fxy);

% Set Seed for random number generator.
rng(SETTINGS.SEED)

% Generate random variables r_{i}
a = -1;
b = 1;
r = a + ((b - a).* rand(m+1,m+1)); 


% Get vector 'eps' - the vector which stores the upper Noise/Signal
% threshold \epsilon_{i} for each cofficient a_{i}
a = el;
b = eu;
eps = a + rand(m + 1, m + 1).*(b - a);

% Calculate the noise vector = a_{i}r_{i}.*\epsilon_{i}
noise_matrix = fxy .* r .* eps;

% Calculate the perturbed coefficients = a_{i} + a_{i}r_{i}\epsilon_{i}
f_noisy = fxy + noise_matrix;




end
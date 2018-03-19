
function [vNoise, fx_noisy] = noise(fx, v, startrng)

% This function adds noise in the componentwise sense to the coefficients 
% of the Bernstein basis polynomial whose coefficients are stored in the
% vector P. The vector p stores the noisy coefficients of the polynomial.

% P        : vector of coefficients of a Bernstein polynomial

% v        : The ratio (noise level)/(signal level)

% startrng : a positive integer for the random number generator. It is
%            included to make sure that different noise samples are 
%            used for each call to this function.

% p        : The noisy polynomial, that is, the polynomial after
%            noise has been added to P

% noisep   : The vector of the noise added to P

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=length(fx);  % the degree of the polynomial is s-1

% Define the vector of noise samples in the range [-1,...,+1].
% Since this function is used for both polynomials, the seed must be a
% constant that takes a different value for the two polynomials, such 
% that a different seed is used for each polynomial.
rand('seed',startrng);
rp=(2*rand(1,s))-ones(1,s);

% Form p, the vector of coefficients of the noisy polynomial.
vNoise = fx.*rp*v;
fx_noisy = fx + vNoise; 


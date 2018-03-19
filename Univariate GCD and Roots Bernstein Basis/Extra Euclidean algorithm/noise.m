
function [p]=noise(P,v,startrng)

% This function adds noise in the componentwise sense to the coefficients 
% of the Bernstein basis polynomial whose coefficients are stored in the
% vector P. The vector p stores the noisy coefficients of the polynomial.

% startrng is a positive integer for the random number generator. It is
% included to make sure that different noise samples are used for
% each call to this function.

% The noise level is v, which is defined as (noise level)/(signal level).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=length(P);  % the degree of the polynomial is s-1

% Define the vector of noise samples in the range [-1,...,+1].
% Since this function is used for both polynomials, the seed must be a
% constant that takes a different value for the two polynomials, such 
% that a different seed is used for each polynomial.
rand('seed',startrng);
rp=(2*rand(1,s))-ones(1,s);

% Form the vector of noisy coefficients.
p=P+(P.*rp*v);


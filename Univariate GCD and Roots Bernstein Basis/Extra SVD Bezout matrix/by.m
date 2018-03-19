
function [z]= by(m,k)

% This function calculates the binomial coefficient mCk = m choose k.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z=factorial(m)/(factorial(m-k)*factorial(k));

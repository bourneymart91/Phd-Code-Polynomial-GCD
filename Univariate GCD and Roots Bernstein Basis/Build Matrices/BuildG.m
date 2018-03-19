function G = BuildG(t)

% Build the matrix G, which is a diagonal matrix whose entries correspond
% to the binomial coefficients of the GCD d(x). 
%
% Inputs.
%
%
% t :   Degree of GCD

% Build G


G = diag(GetBinomials(t));

end

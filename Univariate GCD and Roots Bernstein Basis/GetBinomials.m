function Bi_m = GetBinomials(m)
% Get the set of binomials \binom{m}{i} for i = 0,\dots,m
%
% % Inputs
%
% m : (Int)
%
% % Outputs
%
% Bi_m : (Vector) Set of binomial coefficients \binom{m}{i}

% Initialise vector to store binomials
Bi_m = zeros(m+1,1);

% Get binomial coefficients and store in vector
for i = 0 : 1 : m
    Bi_m(i+1) = nchoosek(m, i);
end


end
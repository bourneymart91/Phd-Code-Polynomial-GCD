function Q1 = BuildQ1(n_k)
% BuildQ1(n_k)
%
% Build the partition of the matrix Q, corresponding to the binomial
% coefficients of v(y). (Note this function works for both v(y) and u(y)
%
% % Inputs
%
% n : degree of polynomial g(y)
%
% t : index of subresultant S_{t} to be built.

% Produce the matrix Q1
Q1 = diag(GetBinomials(n_k));

end
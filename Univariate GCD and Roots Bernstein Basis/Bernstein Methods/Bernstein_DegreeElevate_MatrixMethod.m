function [c] = Bernstein_DegreeElevate_MatrixMethod(fx,r)
% Given the coefficients of the polynomial f(x) in Bernstein basis, and 
% the number of degree elevations r. Perform r fold degree elevation, by
% matrix multiplictation. 
%
% Inputs.
%
% fx : coefficients of polynomial f(x)
%
% r :  number of degree elevations so that output polynomial is of degree 
%      m+r


% Get the degree of polynomial f(x).
m = size(fx,1) - 1;

% Get diagonal matrix D^{-1}.
D = diag(1./GetBinomials(m+r));

% Initialise matrix E
E = zeros(m+r+1,m+1);

% For each column of the matrix E
for j = 0:1:m
    % for each row of the matrix E
    for i = j:1:j+r
        E(i+1,j+1) = nchoosek(r,i-j) * nchoosek(m,j);
    end
end

% Get the coefficients of the degree elevated polynomial.
c = D*E*fx;

end
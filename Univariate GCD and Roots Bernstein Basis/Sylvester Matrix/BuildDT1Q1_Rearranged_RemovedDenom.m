function DT1Q1 = BuildDT1Q1_Rearranged_RemovedDenom(fx,n_k)
% Build D^{-1}_{m+n-k}T_{n-k}(f)Q_{n-k}, the first partition of the 
% Sylvester subresultant matrix S_{k}(f,g). This function removes a common
% denominator of all non-zero entries in the matrix.
%
% % Inputs
%
% fx : Coefficients of the polynomial f(x)
%
% n_k : Degree of polynomial v(x,y) = n - t


% Get Degree of input polynomial
m = GetDegree(fx);

% Initialise the partition of DTQ \in\mathbb{R}^{(m+n-t+1)\times(n-t+1)}.
DT1Q1 = zeros(m+n_k+1,n_k+1);

% for each column k in the matrix D^{-1}T_{n-k}(f)Q_{n-k}.
for j = 0:1:n_k
    % for each coefficient a_{i} in the polynomial f(x)
    for i = j:1:m+j
        DT1Q1(i+1,j+1) = ...
            fx(i-j+1) .*...
            nchoosek(m+n_k-i,m-(i-j)) .* ...
            nchoosek(i,j);
    end
end

% Common Denominator is included in the coefficient matrix.
%DT1Q1 = DT1Q1 ./ nchoosek(m+n_k,n_k);


end
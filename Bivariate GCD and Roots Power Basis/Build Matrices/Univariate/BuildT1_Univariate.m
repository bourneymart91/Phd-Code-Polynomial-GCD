function T1 = BuildT1_Univariate(fx, n_k)
% Build the (n-k)th convolution matrix for polynomial f(x)
%
% % Inputs
%
% fx : (Vector) Coefficient vector of polynomial f(x)
%
% n_k : (Int) Degree of polynomial v(x)
%
% % Outputs
%
% T1 : (Matrix) Convolution matrix of the univariate polynomial f(x)s

% Get degree of f(x)
m = GetDegree_Univariate(fx);

% Initialise convolution matrix
T1 = zeros(m + n_k + 1, n_k+1);

% Each column of T_{n-k}(f(x))
for i = 0:1:n_k
    
    T1(i+1:m+i+1,i+1) = fx;
    
end

end
function [vMax_ai,vMin_ai] = GetMaxMin(fx, n_k)
% Get the maximum and minimum occurence of a_{i} in C_{n-k}(f) for all
% a_{i}, where a_{i} are the m+1 coefficients of the polynomial f(x)
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% n_k : Subscript of T_{n-k}(f(x))
%
% % Outputs
%
% vMax_ai : (Vector) For each coefficient a_{i}, the vector contains the 
%   maximum occurence of a_{i} in the matrix C_{n-k}(f).
% 
% vMin_ai : (Vector) For each coefficient a_{i}, the vector contains the 
%   minimum occurence of a_{i} in the matrix C_{n-k}(f)
%
%
%
%
% Note : The matrix C_{n-k}(f) can take many forms, called variants





% Get the degree of f(x)
m = GetDegree(fx);

% Get absolute values of f(x).
fx = abs(fx);

% Initialise vectors to store maximum and minimum of each a_{i} in
% T_{n-k}(f)
vMax_ai = zeros(m + 1, 1);
vMin_ai = zeros(m + 1, 1);

% Build the matrix C_{n-k}(f)
Cf = BuildSubresultant_Partition_2Polys(fx, n_k);

nRows_Cf = n_k + 1;

% For each coefficient a_{i} of polynomial f(x)
for i = 0:1:m
    
    % Get all entries in C_{n-k}(f) containing a_{i}
    if (nRows_Cf) > 1
        
        vec_ai = diag(Cf, -i);
        
    else
        
        vec_ai = Cf(i + 1, 1);
        
    end
    
    
    % Get max entry of coefficient a_{i}
    vMax_ai(i + 1) = max(abs(vec_ai));
    
    % Get min entry of coefficient a_{i}
    vMin_ai(i + 1) = min(abs(vec_ai));
    
end


end


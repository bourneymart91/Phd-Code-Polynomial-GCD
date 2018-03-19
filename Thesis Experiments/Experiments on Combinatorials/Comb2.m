function nck = Comb2(n,k)

t1 = factorial(n);
t2 = factorial(k);
t3 = factorial(n-k);
nck = t1./(t2*t3);


end


% Complexity:
% n! performs n-1 operations
% k! performs k-1 operations
% n-k! performs n-k-1 operations
%
% +2 operations to obtain nck
%
% Total number of operations 
% (n-1) + (k-1) + (n-k+1) + 2
%   = 2n + 1
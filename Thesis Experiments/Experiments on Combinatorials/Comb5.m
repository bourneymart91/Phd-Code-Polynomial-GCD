function nck = Comb5(n,k)

d = gcd(n,k);
q = k/d;

if (k==0)
    nck = 1;
else
    nck = (Comb5(n-1,k-1)/q)*n/d;
    
end
end

% Computational complexity


% comb5 is called k times 
% each comb5 contains
% Computing the GCD = ?? operations
% Computing q : 1 operation
% two divisions and one product = 3 operations.
%
% Total k( ?? + 4)

% Similar to comb3 but:
% -- has added cost of GCD computation
% ++ more robust for various values of n and k
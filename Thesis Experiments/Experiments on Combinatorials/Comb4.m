function nck = Comb4(n,k)

if k<n-k 
    t1 = 1;
    for i =n:-1:n-k+1
        t1 = t1 * i; 
    end
    t2 = factorial(k);
    nck = t1/t2;
else
    t1 = 1;
    for i = n:-1:k+1
        t1 = t1 * i;
    end
    t2 = factorial(n-k);
    nck = t1/t2;

end

% Computational complexity analysis

% if k < n-k
% for loop contains n-(n-k+1) +1  = k operations (products)
% factorial contains k operations (products)
% computing nck +1 operations

% Total = 2k+1

% if k > n-k
% for loop contains n-(k+1) + 1 = n-k operations
% factorial contains n-k operations
% computing nck +1 operations

% Total = 2(n-k) + 1

% So computational complexity by this method is either 2k+1 or 2(n-k) + 1
% depending on whether k or n-k is smaller.

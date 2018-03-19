function nck = Comb3(n,k)

if (k==0)
    nck = 1;
else
    nck = Comb3(n-1,k-1)*(n/k);
end

end


% Complexity
% Comb3 is called k times - We only exit when k = 0.
% Each comb3 has 2 operations, a product and a division

% Number of operations = k*(2)
% Computational Complexity =  O(k)
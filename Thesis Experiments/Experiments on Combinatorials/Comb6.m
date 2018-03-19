function nck = Comb6(n,k)

t = 1;

if k < n-k
    for i = n:-1:n-k+1
        t = t* i/(n-i+1);
    end
else
    for i = n:-1:k+1
        t = t*i/(n-i+1);
    end    
end

nck = t;

end


% Computational Complexity Analysis
%
% if k < n-k
% each loop is performed n-(n-k+1)+1 = k times
% each loop contains 4 operations
%
% if k >= n-k
% each loop is performed n-(k+1)+1 = n-k times
% each loop contains 4 operations


function C1 = BuildDTQ_log_partition(fx, n_k)


m = GetDegree(fx);

C1 = zeros(m + n_k + 1, n_k + 1);

for i = 0 : 1 : m
    
    for j = 0 : 1 : n_k
        
        
        myValue = GetValue(fx(i+1), m, n_k, i, j);
        
        C1(i + j + 1, j + 1) = myValue;
        
    end
end

end


function myValue = GetValue(ai, m, n_k, i, j)

myValue = GetValue_d(ai, m, n_k, i, j);

end

function myValue = GetValue_a(ai, m, n_k, i , j)

myValue = ai *  nchoosek(m,i) * nchoosek(n_k, j) ./ nchoosek(m + n_k, i + j);

end

function myValue = GetValue_b(ai, m, n_k, i , j)

myValue = ai *  ...
    exp(...
    gammaln(i+j) ...
    + gammaln(m + n_k - i - j) ...
    + gammaln(m) ...
    + gammaln(n_k) ...
    - gammaln(i) ...
    - gammaln(j) ...
    - gammaln(m - i) ...
    - gammaln(n_k - j) ...
    - gammaln(m + n_k));

end

function myValue = GetValue_c(ai, m, n_k, i , j)

myValue = ai *  ...
    10^(...
    log(nchoosek(m, i)) ...
    + log(nchoosek(n_k, j)) ...
    - log(nchoosek(m + n_k, i + j)) ...
    )

end

function myValue = GetValue_d(ai, m, n_k, i , j)

myValue = ai *  ...
    10^(...
    log(factorial(i+j)) ...
    + log(factorial(m + n_k - i - j)) ...
    + log(factorial(m)) ...
    + log(factorial(n_k)) ...
    - log(factorial(i)) ...
    - log(factorial(j)) ...
    - log(factorial(m - i)) ...
    - log(factorial(n_k - j)) ...
    - log(factorial(m + n_k)));

end





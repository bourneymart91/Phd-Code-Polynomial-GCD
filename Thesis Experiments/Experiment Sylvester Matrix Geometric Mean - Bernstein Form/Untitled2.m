
temp_prod = 1;
% Get the product \binom{m+n-k}{i+j}
% For each column in C_{n-k}(f)
for j = 0:1:n_k
    temp_prod = temp_prod .* nchoosek(m+n_k,i+j);
end

test1 = temp_prod;

% Test2
% Get the product \binom{m+n-k+1}{i+j} * \frac{}{}
temp_prod = 1;
for j = 0:1:n_k
   temp_prod = temp_prod.* (nchoosek(m+n_k+1,i+j) * (m+n_k-i-j+1)/(m+n_k+1))  ;
end

test2 = temp_prod;

%test3 
prod3a = 1;
for j=0:1:n_k+1
    prod3a = prod3a * nchoosek(m+n_k+1,i+j)
end

prod3c = 1;
for j = 0:1:n_k
    prod3c = prod3c * ((m+n_k-i-j+1)/(m+n_k+1));
end

prod3b = 1./nchoosek(m+n_k+1,i+n_k+1);

test3 = prod3a * prod3b * prod3c;

test3./test2
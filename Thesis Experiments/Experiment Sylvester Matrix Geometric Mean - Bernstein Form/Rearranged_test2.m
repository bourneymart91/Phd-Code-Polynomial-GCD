function [] = Rearranged_test2(m,n_k)
% Testing geometric means of component parts.
%
% 1. Geometric mean of a_{i} in S_{k} and S_{k-1}
%
% 2. Geometric mean of \binom{i+j}{i} in S_{k} and S_{k-1}
%
% 3. Geometric mean of \binom{m+n-k-i-j}{} in S_{k} and S_{k-1}
%
% 4. Geometric mean of \binom{m+n-k}{m} in S_{k} and S_{k-1}

% Set i to be any value between 0 and m.
i = 1;

% %
% % 
% Test part 2
% 2. Geometric mean of \binom{i+j}{i} in S_{k} and S_{k-1}

% Get geometric mean of first binomial in kth subresultant
temp_prod = 1;
for j = 0:1:n_k
   temp_prod = temp_prod .* nchoosek(i+j,i); 
end

GM_kth = temp_prod.^(1./(n_k+1));


% Get Geometric mean of first binomial in k-1 th subresultant
temp_prod = 1;
for j = 0:1:n_k+1
   temp_prod = temp_prod.* nchoosek(i+j,i);
end

GM_k_minus_1th = temp_prod.^(1./(n_k+2));

% Get geometric mean of first binomial in kth subresultant using the
% geometric mean of the first binomial in the k-1 th subresultant.
test = ( (GM_k_minus_1th.^(n_k+2)) .* (1./nchoosek(i+n_k+1,i)) ) .^(1./(n_k+1))

% %
% %
% % 
% Test Part 3.
% 3. Geometric mean of \binom{m+n-k-i-j}{} in S_{k} and S_{k-1}

% Get geometric mean of kth subresultant
temp_prod = 1;
for j = 0:1:n_k
   temp_prod  = temp_prod .* nchoosek(m+n_k-i-j,m-i);
end
GM_kth = temp_prod.^(1./(n_k+1));

% Get geometric mean of k-1 th subresultant
temp_prod = 1;
for j = 0:1:n_k+1
   temp_prod = temp_prod .* nchoosek(m+n_k-i-j+1,m-i);
end
GM_k_minus_1th = temp_prod.^(1./(n_k+2));

% Test obtaining GM(Sk) from GM(Sk-1)

% Get Geometric mean of k-th from k-1th and test
test = (GM_k_minus_1th.^(n_k+2)) .* (1./(nchoosek(m+n_k-i+1,m-i)));
test = test.^(1./(n_k+1));


% %
% %
% %
% Test Part 4
% 4. Geometric mean of \binom{m+n-k}{m} in S_{k} and S_{k-1}

% Get geometric mean of the denominator in the kth subresultant 
temp_prod = 1;
for j = 0:1:n_k
    temp_prod = temp_prod .* (1./ nchoosek(m+n_k,m));
end
GM_kth = temp_prod.^(1./(n_k+1))

% Get geometric mean of the denominator in the k-1 th subresultant
temp_prod = 1;
for j = 0:1:n_k+1
   temp_prod = temp_prod .* (1./nchoosek(m+n_k+1,m));
end
GM_k_minus_1th = temp_prod.^(1./ (n_k+2))

% Get Geometric mean of the denominator in the kth subresultant from
% geometric mean of the denominators in the  k-1th subresultant.
frac = ((m+n_k+1) ./ (n_k+1))^(n_k+2)
GM_test = (GM_k_minus_1th.^(n_k+2) * nchoosek(m+n_k,m) * frac).^(1./(n_k+1))

end
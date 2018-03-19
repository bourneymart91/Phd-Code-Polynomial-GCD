function [] = Standard_test1(m,n_k)
% Testing the relationship between product of individual elements in
% C_{k}(f) and C_{k+1}(f).
%
% 1. Product of a_{i} in S_{k} and product of a_{i} in S_{k-1}
%
% 2. Product of \binom{m}{i} in S_{k} and S_{k-1}
%
% 3. Product of \binom{n-k}{j} in S_{k} and S_{k-1}
%
% 4. Product of \frac{1} {\binom{m+n-k}{i-j}} in S_{k} and S_{k-1}

%test1 

% Pick i between 0 and m
i = 1;

% % 
% % 
% %
% %
% Test One 

% %
% %
% %
% %
% Test Two

temp_prod = 1;
for j = 0:1:n_k
    temp_prod = temp_prod * nchoosek(m,i);
end
test2a = temp_prod

temp_prod = 1;
for j = 0:1:n_k+1
    temp_prod = temp_prod * nchoosek(m,i);
end

test2b = temp_prod * (1./ (nchoosek(m,i)));

display(test2a - test2b)

% %
% % 
% %
% %
% Test 3

% Get product of second binomial in kth subresultant

temp_prod = 1;
for j = 0:1:n_k
    temp_prod = temp_prod * nchoosek(n_k,j);
end
test3a = temp_prod;

% Now test derivation from product in terms of k-1 th subresultant
temp_prod = 1;
for j = 0:1:n_k+1
    temp_prod = temp_prod * nchoosek(n_k+1,j);
end
b_part_a = temp_prod;

b_part_b = 1./(nchoosek(n_k+1,n_k+1));

temp_prod = 1;
for j = 0:1:n_k
   temp_prod = temp_prod .* ((n_k-j+1) ./ (n_k+1));
end
b_part_c = temp_prod;

test3b = b_part_a * b_part_b * b_part_c;

display(test3a - test3b)

% % 
% %
% %
% %
% Test 4
temp_prod = 1;
for j = 0:1:n_k
    temp_prod = temp_prod * nchoosek(m+n_k,i+j);
end
test4a = temp_prod;

temp_prod = 1;
for j = 0:1:n_k
   temp_prod = temp_prod .* (nchoosek(m+n_k+1,i+j) * ((m+n_k-i-j+1)./(m+n_k+1)));
end
test4b = temp_prod

temp_prod = 1;
for j = 0:1:n_k+1
    temp_prod = temp_prod * nchoosek(m+n_k+1,i+j)
end
b_part_a = temp_prod;

b_part_b = 1./nchoosek(m+n_k+1,i+n_k+1);

temp_prod = 1;
for j = 0:1:n_k
    temp_prod = temp_prod * ((m+n_k-i-j+1)./(m+n_k+1));
end
b_part_c = temp_prod;

test4b = b_part_a * b_part_b * b_part_c;

display(test4a-test4b)


end
function [] = Standard_test2(m,n_k)


i = 1;

% %
% % 
% %
% Test Part One


% % 
% % 
% % 
% Test Part Two

% Get Geometric mean of the first binomial in kth subresultant
prod = 1;
for j = 0:1:n_k
    prod = prod * nchoosek(m,i);
end

test2 = prod.^(1./(n-k+1));

% Get Geometric mean 
prod = 1;
for j = 0:1:n_k+1
    prod = prod * nchoosek(n_k+1,j);
end

test_2a = prod

test_2b = 1./(nchoosek(n_k+1,n_k+1);

prod = 1;
for j = 0:1:n_k
    prod = prod .* (n_k-j+1)./(n_k+1);
end

test_2c = prod;



end

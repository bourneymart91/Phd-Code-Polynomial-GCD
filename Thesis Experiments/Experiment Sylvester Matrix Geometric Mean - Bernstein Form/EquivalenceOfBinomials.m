
% prod1 : The product of the first binomial coefficient in the numerator of
% the entries of C(f)

% prod2 : The product of the second binomial coefficient in the numerator
% of the entries of C(f)

n_k = 10;
m = 5;

prod = 1;
for i = 0:1:m
    for j = 0:1:n_k
        prod = prod * nchoosek(i+j,i);
    end
end

prod2 = 1;
for i = 0:1:m
    for j = 0:1:n_k
        prod2 = prod2 * nchoosek(m+n_k-(i+j),m-i) ;
    end
end

display(prod)
display(prod2)
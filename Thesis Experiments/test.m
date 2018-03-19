
% Experiment on the sum of the diagonals of the convolution matrix
% C_{n-k}(f(x,y))

% Set degree of polynomial f(x,y)
m = 5;

% Set degree of polynomial v(x,y)
n_k = 22;


% Initialise vector to store the sum of each a_{i}
sum = zeros(m+1,1);
prod = ones(m+1,1);

for i = 0:1:m
    for j = 0:1:n_k

        sum(i+1) = sum(i+1) + (nchoosek(i+j,i) * nchoosek(m+n_k-i-j,m-i));
        prod(i+1) = prod(i+1) * ((nchoosek(i+j,i) * nchoosek(m+n_k-i-j,m-i)));
    end

end



display(sum);
display(prod);

function Q1 = BuildQ1(m)
% BuildQ1(m)
%
% Given the coefficients of polynomial f(x,y), build the diagonal matrix Q
% of trinomial coefficients. 
% Used in convolution the product f*v is given by 
% D^{-1}T_{n-k}(f)Q_{n-k} * v = h.
% Also used in construction of the kth Sylvester subresultant matrix. The
% matrix S_{k}(f,g) = D^{-1}T(f,g)blkdiag(Q_{n-k},Q_{m-k}) 
%
% Inputs
%
% m : Degree of polynomial f(x,y)
%
% Outputs.
%
% Q : The diagonal matrix Q1

% Initialise a temporary matrix to store the trinomial coefficients.
temp_mat = zeros(m+1,m+1);

for i = 0:1:m
    for j = 0:1:m-i
    
        temp_mat(i+1,j+1) = Trinomial(m,i,j);
    
    end
end

% Get the entries as a vector
temp_vec = GetAsVector(temp_mat);

% Get number of nonzeros in the temp vector
nNonZeros = nchoosek(m+2,2);

% Remove the zeros from the temp vector
temp_vec = temp_vec(1:nNonZeros);

% Build the matrix Q_{m}
Q1 = diag(temp_vec);


end
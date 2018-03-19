function D = BuildD_2Polys(m, n_k)
% BuildD(m, n_t)
%
% Build the diagonal matrix D^{-1} for the convolution of two polynomials
% f(x,y) and v(x,y) of degrees m and n-t.
%
% Inputs
%
% m : (Int) Degree of polynomial f(x,y)
%
% n_k : (Int) Degree of polynomial v_{k}(x,y)

% Initialise a temporary matrix to store the (m+n-t+1) x (m+n-t+1) diagonal
% entries of the matrix D.
temp_mat = zeros(m+n_k+1, m+n_k+1);

% Fill entries of temp matrix with trinomial coefficients
for i = 0:1:m+n_k
    
    for j = 0:1:m+n_k-i    
        
        temp_mat(i+1,j+1) = 1./Trinomial(m + n_k, i, j);
        
    end
    
end

% Read entries of the temp matrix into a vector
vect = GetAsVector(temp_mat);

% Only want the nchoosek(m+n+2,2) non-zero terms
nNonZeroTerms = nchoosek(m+n_k+2, 2);

% Remove the zero terms.
vect = vect(1 : nNonZeroTerms);

% Get diagonal matrix D^{-1}
D = diag(vect);

end

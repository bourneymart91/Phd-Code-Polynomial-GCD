function P1 = BuildP1(m, n, k, idx_col)
% Build the matrix G1, which forms a partition of the matrix P, and is used
% in the BuildP() function. This function forms part of the STLN package.
%
% % Input
%
% m : (Int) Degree of polynomail f(x,y)
% 
% n : (Int) Degree of polynomial g(x,y)
%
% k : (Int) Degree of polynomial d(x,y) and index of Sylvester subresultant
% S_{k}(f,g)
%
% idx_col : (Int) index of column removed from S_{k}(f,g)
%
% % Outputs
%
% P1 : (Matrix) P1


% Get the number of columns in the first partition of the Sylvester matrix
nCoefficients_vxy = nchoosek(n-k+2,2);
nColumns_Tf = nCoefficients_vxy;

% Get the number of zeros in v(x,y)
nZeros_vxy = nchoosek(n-k+1,2);


vec = (1:1:nColumns_Tf)';
vec = [vec ; zeros(nZeros_vxy,1)];

mat = GetAsMatrix(vec, n-k, n-k);

[row,col] = find(mat==idx_col);

row_index = row - 1;
col_index = col - 1;

j1 = row_index;
j2 = col_index;

trinom = Trinomial(n-k, j1, j2);

% Build a matrix of ones corresponding to f(x,y)
nCoefficients_fxy = nchoosek(m+2,2);
nZeros_fxy = nchoosek(m+1,2);

fxy = [ones(nCoefficients_fxy,1); zeros(nZeros_fxy,1)];
fxy = GetAsMatrix(fxy,m,m);

% Get matrix of the product of f(x,y) and v(x,y)
prod = zeros(m+n-k+1,m+n-k+1);

prod(j1+1:j1+m+1,j2+1:j2+m+1) = fxy;

vec = GetAsVector(prod);

% remove zeros
nCoefficients_prod = nchoosek(m+n-k+2,2);
nZeros_prod = nchoosek(m+n-k+1,2);

vec = vec(1:nCoefficients_prod);

A = diag(vec);

% Remove the zero columns
A(:, find(sum(abs(A)) == 0)) = [] ;
P1 = A .* trinom;
end
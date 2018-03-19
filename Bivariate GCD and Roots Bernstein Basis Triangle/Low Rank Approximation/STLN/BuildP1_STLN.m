function P1 = BuildP1_STLN(m,n,k,idx_col)
% Build the matrix G1, which forms a partition of the matrix P, and is used
% in the BuildP() function. This function forms part of the STLN package.
%
% % Inputs
%
% m : (Int) Total degree of polynomial f(x,y)
%
% n : (Int) Total degree of polynomial g(x,y)
%
% k : (Int) Total degree of polynomial d(x,y) and index of Sylvester subresultant
%       S_{k}(f,g) whose low rank approximation is computed by STLN.
%
% % Outputs
%
% P1

% Get the number of columns in the first partition of the Sylvester matrix
nColumns_Tf = nchoosek(n-k+2,2);

% Get the number of zeros
try
    nZeros = nchoosek(n-k+1,2);
catch
    nZeros = 0;
end

vec = (1 : 1 : nColumns_Tf)';
vec = [vec ; zeros(nZeros,1)];

mat = GetAsMatrix(vec, n-k, n-k);

[row,col] = find(mat==idx_col);

row_index = row - 1;
col_index = col - 1;

j1 = row_index;
j2 = col_index;

trinom = Trinomial(n-k,j1,j2);

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
nCoefficients_fv = nchoosek(m+n-k+2,2);
nZeros_prod = nchoosek(m+n-k+1,2);

vec = vec(1:nCoefficients_fv);

A = diag(vec);

% Remove the zero columns
A(:, find(sum(abs(A)) == 0)) = [] ;
P1 = A .* trinom;
end
function P = BuildP1_TotalDegree_SNTLN(m,n,k,idx_col)
% Build the matrix P, used in SNTLN function. P is a matrix which is
% obtained from the decomposition of a column vector c_{t} into a matrix
% vector product P_{t} [f;g]
% c_{t} is the column of the Sylvester matrix S_{t}(f,g).
%
% %   Inputs
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% k : (Int) Total degree of d(x,y)
%
% th1 : (Float) Optimal value of theta_{1}
%
% th2 : (Float) Optimal value of theta_{2}
%
% idx_col : (Int) Optimal column for removal from S(f,g)


% Get the number of coefficients in f(x,y)
nNonZeros_fxy = nchoosek(m+2,2);
% Get the number of zeros in the matrix f(x,y)
nZeros_fxy = nchoosek(m+1,2);


% Get a matrix the same size as f(x,y) where the coefficients are replaced
% by ones.
vec = [...
    ones(nNonZeros_fxy,1);...
    zeros(nZeros_fxy,1);...
    ];
mat = GetAsMatrix_Version1(vec,m,m);


% Produce a zero matrix to fill the space
padd_mat = zeros(m+n-k+1, m+n-k+1);


% % from the index of the optimal column, Get the number of multiplications
% with respec to x and number with respect to y.
[i,j] = GivenCol_GetIndex_Relative(n-k,n-k,idx_col);


ihat = i+1;
jhat = j+1;

nRows_f = m+1;
nCols_f = m+1;

% inser the theta matrix into the zero matrix
padd_mat(ihat:i+nRows_f, jhat:j+nCols_f) = mat;

% Get the padded matrix as a vector
vec_padd_mat = GetAsVector_Version1(padd_mat);
nCoeff_fv = nchoosek(m+n-k+2,2);
vec_padd_mat = vec_padd_mat(1:nCoeff_fv);


% Diagonalise the vector.
diag_mat_vec_padd_mat = diag(vec_padd_mat);


% Remove the zero columns
diag_mat_vec_padd_mat( :, ~any(diag_mat_vec_padd_mat,1) ) = [];  %columns


P = diag_mat_vec_padd_mat;


end
function P1 = BuildP1_BothDegree_SNTLN(m,m1,m2,n,n1,n2,k,k1,k2,idx_col)
% BuildP1_BothDegree_SNTLN(m,m1,m2,n,n1,n2,k,k1,k2,idx_col)
%
% Build the matrix P, used in SNTLN function. P is a matrix which is
% obtained from the decomposition of a column vector c_{t} into a matrix
% vector product P_{t} [f;g]
% c_{t} is the column of the Sylvester matrix S_{t}(f,g).
%
% % Inputs
%
% m : (Int) Total degree of f(x,y)
%
% m1 : (Int) Degree of polynomial f(x,y) with respect to x
%
% m2 : (Int) Degree of polynomial f(x,y) with respect to y
%
% n : (Int) Total degree of g(x,y)
%
% n1 : (Int) Degree of polynomial g(x,y) with respect to x
%
% n2 : (Int) Degree of polynomial g(x,y) with respect to y
%
% k : (Int) Total degree of d(x,y)
%
% k1 : (Int) Degree of GCD d(x,y) with respect to x
%
% k2 : (Int) Degree of GCD d(x,y) with respect to y
%
% th1 : (Float) Optimal value of theta_{1}
%
% th2 : (Float) Optimal value of theta_{2}
%
% idx_col : (Int) Optimal column for removal from S(f,g)
%
% % Outputs
%
% P1 : (Matrix) 

% Get number of nonzeros of f(x,y).
nNonZeros_fxy = GetNumNonZeros(m1,m2,m);
nZeros_fxy_relative = ((m1+1) * (m2+1)) - nNonZeros_fxy;

% Build a matrix the same size of f(x,y) with the coefficients replaced by
% ones.
vec = [...
    ones(nNonZeros_fxy,1);
    zeros(nZeros_fxy_relative,1)
    ];

mat = GetAsMatrix_Version1(vec,m1,m2);


% Produce a zero matrix.
padd_mat = zeros(m1+n1-k1+1, m2+n2-k2+1);

% % From the index of the optimal column, Get the number of multiplications
% with respec to x and number with respect to y.
[i,j] = GivenCol_GetIndex_Relative(n1-k1, n2-k2, idx_col);

ihat = i+1;
jhat = j+1;

nRows_f = m1+1;
nCols_f = m2+1;

% Insert the theta matrix into the zero matrix
padd_mat(ihat:i+nRows_f, jhat:j+nCols_f) = mat;

% Get the padded matrix as a vector
vec_padd_mat = GetAsVector_Version1(padd_mat);

% Get the number of coefficients in the product fv
nCoeff_fv = GetNumNonZeros(m1+n1-k1,m2+n2-k2,m+n-k);

% Diagonalise the vector.
diag_mat_vec_padd_mat = diag(vec_padd_mat);

% Remove the zero columns
diag_mat_vec_padd_mat( :, ~any(diag_mat_vec_padd_mat,1) ) = [];  %columns

P1 = diag_mat_vec_padd_mat;
P1 = P1(1:nCoeff_fv,:);

end
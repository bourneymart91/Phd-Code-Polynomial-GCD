function Pk = BuildP_BothDegree_SNTLN(m,m1,m2,n,n1,n2,k,k1,k2,alpha,th1,th2,idx_col)
% Calculate the matrix P where P is the matrix such that a column of the
% Sylvester subresultant matrix S_{k,k1,k2} can be written as a product of
% P and the column vector of coefficients of f and g.
%
% Inputs
%
% m : (Int) Total degree of polynomial f(x,y)
%
% m1 : (Int) Degree of f(x,y) with respect to x
%
% m2 : (Int) Degree of f(x,y) with respect to y
%
% n : (Int) Total degree of polynomial g(x,y)
%
% n1 : (Int) Degree of g(x,y) with respect to x
%
% n2 : (Int) Degree of g(x,y) with respect to y
%
% k : (Int) Total degree of d(x,y)
%
% k1 : (Int) Degree of d(x,y) with respect to x
%
% k2 : (Int) Degree of d(x,y) with respect to y
%
% alpha : (Float) Optimal value of \alpha
%
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}
%
% idx_col : (Int) Index of column removed from S_{k_{1},k_{2}}
%
%
% % Outputs
%
% P : (Matrix) P

% Get the number of coefficients in polynomial f
% % nCoeff_f = (m1+1).*(m2+1);
nNonZeros_fxy = GetNumNonZeros(m1, m2, m);
nZeros_fxy = (m1+1)*(m2+1) - nNonZeros_fxy;

% Get the number of coefficients in polynomial g(x,y)
% % nCoeff_g = (n1+1).*(n2+1);
nNonZeros_gxy = GetNumNonZeros(n1, n2, n);
nZeros_gxy = (n1+1)*(n2+1) - nNonZeros_gxy;

% Get number of Rows in Sylvester matrix S_{k,k1,k2}
nRows_Skk1k2 = GetNumNonZeros(m1+n1-k1,m2+n2-k2,m+n-k);


% Get number of columns in first partition of the Sylvester subresultant
% matrix.
nColumns_T1 = GetNumNonZeros(n1-k1,n2-k2,n-k);

if idx_col <= nColumns_T1
    % Optimal column in first partition
    
    % % Build the matrix P
    
    % Build the matrix P1
    P1 = BuildP1_BothDegree_SNTLN(m, m1, m2, n, n1, n2, k, k1, k2, idx_col);
    
    % Build the matrix P2
    P2 = zeros(nRows_Skk1k2,nNonZeros_gxy);
    
else
    % Optimal column in second partition
    
    % Build the matrix P1
    P1 = zeros(nRows_Skk1k2,nNonZeros_fxy);
    
    % Build the matrix P2
    % Get the position of the optimal column with respect to T(g)
    idx_col_rel = idx_col - nColumns_T1;
    
    P2 = BuildP1_BothDegree_SNTLN(n, n1, n2, m, m1, m2, k, k1, k2, idx_col_rel);
    
    
    
end

f_vec = [ones(nNonZeros_fxy,1); zeros(nZeros_fxy,1)];
f_mat = GetAsMatrix(f_vec,m1,m2);

g_vec = [ones(nNonZeros_gxy,1); zeros(nZeros_gxy,1)];
g_mat = GetAsMatrix(g_vec,n1,n2);


% Get a vector of theta_{1}theta_{2} corresponding to entries of the
% polynomial f(x,y)
th_f = GetWithThetas(f_mat,th1,th2);
%th_f = GetWithThetas(f_mat,1,1);
th_f = GetAsVector(th_f);
th_f = th_f(1:nNonZeros_fxy);
th_f = diag(th_f);

% Get a vector of \theta_{1}\theta_{2} corresponding to entires of the
% polynomial g(x,y)
th_g = GetWithThetas(g_mat,th1,th2);
%th_g = GetWithThetas(g_mat,1,1);
th_g = GetAsVector(th_g);
th_g = th_g(1:nNonZeros_gxy);
th_g = diag(th_g);

% Build the matrix P.
Pk = [P1*th_f alpha.*P2*th_g];

end

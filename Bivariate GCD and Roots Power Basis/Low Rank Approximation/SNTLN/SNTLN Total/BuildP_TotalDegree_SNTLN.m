function P = BuildP_TotalDegree_SNTLN(m, n, k, alpha, th1, th2, idx_col)
% Calculate the matrix P where P is the matrix such that a column of the
% Sylvester subresultant matrix S_{k,k1,k2} can be written as a product of
% P and the column vector of coefficients of f and g.
%
% Inputs
%
% m : Total degree of polynomial f(x,y)
%
% n : Total degree of polynomial g(x,y)
%
% k : Total degree of polynomial d(x,y)
%
% alpha : Optimal value of \alpha
%
% th1 : Optimal value of \theta_{1}
%
% th2 : Optimal value of \theta_{2}
%
% idx_col : Index of column removed from S_{k_{1},k_{2}}
%
% % Outputs
%
% P : Matrix P

% Get the number of coefficients in polynomial f
nNonZeros_fxy = nchoosek(m+2,2);

% Get the number of coefficients in polynomial g(x,y)
nNonZeros_gxy = nchoosek(n+2,2);

% Get number of Rows in Sylvester matrix S_{k,k1,k2}
nRows_Skk1k2 = nchoosek(m+n-k+2,2);


% Get number of columns in first partition of the Sylvester subresultant
% matrix.
nCols_T1 = nchoosek(n-k+2,2);

if idx_col <= nCols_T1 % Optimal column in first partition
    
       
    % Build the matrix P1
    P1 = BuildP1_TotalDegree_SNTLN(m,n,k,idx_col);
    
    % Build the matrix P2
    P2 = zeros(nRows_Skk1k2,nNonZeros_gxy);
    
    
else
    % Optimal column in second partition
    
    % Build the matrix P1
    P1 = zeros(nRows_Skk1k2,nNonZeros_fxy);
    
    % Build the matrix P2
    % Get the position of the optimal column with respect to T(g)
    opt_col_rel = idx_col - nCols_T1;
    
    P2 = BuildP1_TotalDegree_SNTLN(n,m,k,opt_col_rel);
    
    
    
end

% Get thetas corresponding to the coefficients of f(\omega_{1},\omega_{2})
pre_th1s = diag(th1.^(0:1:m));
post_th2s = diag(th2.^(0:1:m));
fww_thetas_mat = pre_th1s * ones(m+1,m+1) * post_th2s;
vec = GetAsVector_Version1(fww_thetas_mat);
vec = vec(1:nNonZeros_fxy);
th_f = diag(vec);

% Get thetas corresponding to the coefficients of f(\omega_{1},\omega_{2})
pre_th1s = diag(th1.^(0:1:n));
post_th2s = diag(th2.^(0:1:n));
gww_thetas_mat = pre_th1s * ones(n+1,n+1) * post_th2s;
vec = GetAsVector_Version1(gww_thetas_mat);
vec = vec(1:nNonZeros_gxy);
th_g = diag(vec);

% Build the matrix P.
P = [P1 * th_f alpha.*P2*th_g];

end
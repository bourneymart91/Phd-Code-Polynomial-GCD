function P = BuildP_RelativeDegree_SNTLN(m1,m2,n1,n2,k1,k2,alpha,th1,th2,idx_col)
% Calculate the matrix DP where P is the matrix such that c = P[f;g]
%
% Inputs
%
% m1 : Degree of f(x,y) with respect to x
%
% m2 : Degree of f(x,y) with respect to y
%
% n1 : Degree of g(x,y) with respect to x
%
% n2 : Degree of g(x,y) with respect to y
%
% k1 : Degree of d(x,y) with respect to x
%
% k2 : Degree of d(x,y) with respect to y
%
% alpha : Optimal value of \alpha
%
% th1 : Optimal value of \theta_{1}
%
% th2 : Optimal value of \theta_{2}
%
% idx_col : Index of column removed from S_{k_{1},k_{2}}
%
% % Outputs.
%
% P : Matrix P

% Get the number of coefficients in polynomial f(x,y)
nCoeffs_fxy = (m1+1).*(m2+1);

% Get the number of coefficients in polynomial g(x,y)
nCoeffs_gxy = (n1+1).*(n2+1);

% Get the number of cols in the first partition of the Sylvester
% subresultant matrix S_{k1,k2}
nCols_T1 = (n1-k1+1) * (n2-k2+1);

nRowsSylvester = (m1+n1-k1+1)*(m2+n2-k2+1);

if idx_col <= nCols_T1 % Optimal column in first partition
    
    
    % % Build the matrix P
    
    % Build the matrix P1
    P1 = BuildP1_RelativeDegree_SNTLN(m1,m2,n1,n2,k1,k2,idx_col);
    
    % Build the matrix P2
    P2 = zeros(nRowsSylvester,nCoeffs_gxy);
    
      
else % Optimal column in second partition
    
    
    % Build the matrix P1
    P1 = zeros(nRowsSylvester,nCoeffs_fxy);
    
    
    % Get the position of the optimal column with respect to T(g)
    opt_col_rel = idx_col - nCols_T1;
    
    % Build the matrix P2
    P2 = BuildP1_RelativeDegree_SNTLN(n1,n2,m1,m2,k1,k2,opt_col_rel);
    
    
    
end

% Get vector of thetas corresponding to the coefficients of f(\omega_{1},\omega_{2})
pre_th1s = diag(th1.^(0:1:m1));
post_th2s = diag(th2.^(0:1:m2));
fww_thetas_mat = ones(m1+1,m2+1);
fww_thetas_mat = pre_th1s * fww_thetas_mat * post_th2s;
th_f = GetAsVector(fww_thetas_mat);
th_f = diag(th_f);

% Get vector of thetas corresponding to the coefficients of g(\omega_{1},\omega_{2})
pre_th1s = diag(th1.^(0:1:n1));
post_th2s = diag(th2.^(0:1:n2));
gww_thetas_mat = ones(n1+1,n2+1);
gww_thetas_mat = pre_th1s * gww_thetas_mat * post_th2s;
th_g = GetAsVector(gww_thetas_mat);
th_g = diag(th_g);


% Build the matrix P.
P = [P1*th_f alpha.*P2*th_g];

end

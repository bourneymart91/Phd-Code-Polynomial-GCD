
function P = BuildP_SNTLN(m, n, k, alpha, theta, idx_col)
% The matrix P is such that a column c_{t} can be expressed as P_{t}[f;g]
% Given a column of the Sylvester matrix S_{t}(f,g), obtain a decomposition
% so that it is expressed as a matrix vector product where the vector
% contains only coefficients of f and g.
%
% % Inputs
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% k : Degree of polynomial d(x)
%
% alpha : Optimal value of \alpha
%
% theta : Optimal value of \theta
% 
% idx_col : Index of optimal column S_{k}(f,g)
%
% % Outputs
%
% P : Matrix P

nCols_Tf = n-k+1;

if idx_col <= nCols_Tf % column is in the 1st partiton of Sk(f,g)
    
    % Build the partiton P1
    P1 = BuildP1(m,n,k,idx_col);
    
    % Build P2
    P2 = zeros(m+n-k+1,n+1);
    
    
else % Column is in the second partition of Sk(f,g)
    
    % Get index of column relative to T_{m-k}(g) (the second partition)
    opt_col_rel = idx_col - (n-k+1);
    
    % Build P1
    P1 = zeros(m+n-k+1,m+1);
    
    % Build P2
    P2 = BuildP1(n,m,k,opt_col_rel);
    
end

% Get thetas corresponding to polynomial f(x)
th_f = diag(GetWithThetas(ones(m+1,1),theta));

% Get thetas corresponding to polynomial g(x)
th_g = diag(GetWithThetas(ones(n+1,1),theta));

% Build the matrix P
P = [P1*th_f alpha.*P2*th_g];


end


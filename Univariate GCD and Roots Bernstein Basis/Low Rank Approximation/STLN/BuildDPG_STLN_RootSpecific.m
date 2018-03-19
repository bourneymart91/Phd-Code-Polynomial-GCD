function DPQ = BuildDPG_STLN_RootSpecific(m, n, k, idx_col, alpha, theta)
% BuildDPQ(m,n,k,idx_col)
%
% Build the matrix DP. Build the matrix DP such that DP * [f;g] gives the
% column of the Sylvester subresultant matrix matrix whose index is given
% by idx_col.
%
%
%
% Inputs
%
% m : (Int) Degree of polynomial f(x)
%
% n : (Int) Degree of polynomial g(x)
%
% theta : (Float) Optimal value of \theta
%
% idx_col : (Int) Index of column c_{k} removed from S_{k}(f,g)
%
% k : (Int) Degree of GCD d(x)
%
% Outputs.
%
% DPQ : (Matrix) DPQ

% Get the number of columns in T_{n-k}(f)
nColumns_Tf = m - k ;

% Build the matrix D^{-1}_{m+n-k}
D = BuildD_2Polys(m, m - k - 1);
    
% Build the matrices P_{1} and P_{2}

if idx_col <= nColumns_Tf % Column is in first partition T_{n-k}(f) of S_{k}

    fprintf('First Partition \n')
    
    j = idx_col - 1;

    
    th_m = diag(GetWithThetas(ones(m + 1, 1), theta));
    
    
    P1 = ...
        [
            zeros(j, m + 1);
            th_m * eye(m + 1, m + 1);
            zeros(m - k - j - 1, m + 1);
        ];
    
    Qm = BuildQ1(m);
    
    DPQ = nchoosek(m - k - 1, j) * D * P1  *Qm;
    
else  %  The column is from the second partiton of the Sylvester matrix poly g
    
    fprintf('Second Partition \n')
    
    % Get index of column relative to second partition T_{m-k}(g)
    idx_col = idx_col - (n - k + 1);
    
    j = idx_col - 1;
    
    Q = BuildQ1(m-1);
    
    T = zeros(m , m + 1);
    for i = 1 : 1 : m
        T(i,i:i+1) = [-1 1];
    end
    
    G = Q * T;
    
    
    th_m = diag(GetWithThetas(ones(m,1),theta));
    
    P1 = ...
        [
            zeros(j,m + 1);
            th_m * G 
            zeros(m - k - j, m + 1);
        ];
    
    DPQ = alpha.* m * nchoosek(m - k , j) * D * P1 ;
    
   
end


end

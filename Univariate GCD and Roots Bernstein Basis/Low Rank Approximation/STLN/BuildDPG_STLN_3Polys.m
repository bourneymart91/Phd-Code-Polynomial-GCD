function DPQ = BuildDPG_STLN_3Polys(m, n, o, k, idx_col)
% BuildDPQ(m,n,k,idx_col)
%
% Build the matrix DP. Build the matrix DP such that DP * [f;g] gives the
% column of the Sylvester subresultant matrix matrix whose index is given
% by idx_col.
%
% Used in SNTLN.m
% Used in STLN.m
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
nColumns_FirstPart = n - k + 1;
nColumns_SecondPart = o - k + 1;
nColumns_ThirdPart = m - k + 1;


% Build the matrix D^{-1}_{m+n-k}
D = BuildD_3Polys(m, n-k, o-k);

% Build the matrices P_{1} and P_{2}

P1 = zeros(m + n - k + 1, m + 1);
P2 = zeros(m + n - k + 1, n + 1);
P3 = zeros(m + n - k + 1, o + 1);
P4 = zeros(m + o - k + 1, m + 1);
P5 = zeros(m + o - k + 1, n + 1);
P6 = zeros(m + o - k + 1, o + 1);


if idx_col <= nColumns_FirstPart % Column is in first partition T_{n-k}(f) of S_{k}
    
    % Build the matrix P
    P1 = BuildP1(m, n - k, idx_col) .* nchoosek(n - k, idx_col - 1);
    
    
    
    
    
elseif idx_col <= nColumns_FirstPart + nColumns_SecondPart  %  The column is from the second partiton of the Sylvester matrix poly g
    
    % Get index of column relative to second partition T_{m-k}(g)
    
    idx_col = idx_col - nColumns_FirstPart;

    % Build the matrix P_{2}
    P4 =  BuildP1(m, o - k, idx_col) .* nchoosek(o - k, idx_col - 1);
    
    
else % Third partition
    
    
    idx_col = idx_col - nColumns_FirstPart - nColumns_SecondPart;
    
    P2 = BuildP1(n, m - k, idx_col) .* nchoosek(m - k, idx_col - 1);
    P6 = BuildP1(o, m - k, idx_col) .* nchoosek(m - k, idx_col - 1);
    
end

% Build the matrices Q
Q1 = BuildQ1(m);
Q2 = BuildQ1(n);
Q3 = BuildQ1(o);
Q = blkdiag(Q1, Q2, Q3);

% Build DPG
DPQ = D*[P1 P2 P3 ; P4 P5 P6] * Q ;

end

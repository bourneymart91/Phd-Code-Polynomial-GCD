function DPG = BuildDPG_SNTLN(m,n,k,alpha,th1,th2,idx_col)
% Build the matrix P_{t} such that any column of the Sylvester matrix
% S_{t}(f,g) is expressed as P_{t}*[f;g]
%
% % Inputs
%
% m : (Int) Degree of polynomial f(x,y)
%
% n : (Int) Degree of polynomial g(x,y)
%
% k : (Int) Index of Sylvester subresultant matrix S_{k}(f,g)
%
% alpha : (Float) Optimal value of \alpha
%
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}
%
% idx_col : (Int) Index of column c_{k} to be removed from S_{k}, and formed by 
% the matrix vector product P*[f;g]
%
% % Outputs 
%
% DPG : (Matrix) DPG such that DPG*[f;g] = c_{k}

% Get the number of columns in the matrix C_{n-k}(f)
nColumns_Tf = nchoosek(n-k+2,2);
nColumns_Tg = nchoosek(m-k+2,2);

% Get the number of coefficients in matrix of f(x,y)
nCoefficients_fxy = nchoosek(m+2,2);
nZeros_fxy = nchoosek(m+1,2);

% Get the number of coefficients in matrix of g(x,y)
nCoefficients_gxy = nchoosek(n+2,2);
nZeros_gxy = nchoosek(n+1,2);

% Build the matrix D^{-1}
D = BuildD_2Polys(m, n-k);

% Build The middle 
if idx_col <= nColumns_Tf  % Column is in first partition of S_{k}
    
    % Build G
    G1 = BuildP1(m, n, k, idx_col);
    G2 = zeros(nchoosek(m+n-k+2,2), nchoosek(n+2,2));
    
    
else    % Opt col is in second partition
    
    % Get the index of the column with respect to the second partition of
    % the Sylvester matrix.
    idx_col = idx_col - nColumns_Tf;
    
    % Build G
    G1 = zeros(nchoosek(m+n-k+2,2), nchoosek(m+2,2));
    G2 = BuildP1(n, m, k, idx_col);
    
end

% Build the matrix Q = [Q1 0 ; 0 Q2]
Q1 = BuildQ1(m);
Q2 = BuildQ1(n);

% Get the vector of \theta_{1}\theta_{2} corresponding to coefficients of
% polynomial f(x,y)
f_vec = [ones(nCoefficients_fxy,1);zeros(nZeros_fxy,1)];
mat_fxy = GetAsMatrix(f_vec,m,m);
th_fxy = GetAsVector(GetWithThetas(mat_fxy,m,th1,th2));
th_fxy = th_fxy(1:nCoefficients_fxy);
th_fxy = diag(th_fxy);

% Get the vector of \theta_{1}\theta_{2} corresponding to coefficients of
% polynomial g(x,y)
g_vec = [ones(nCoefficients_gxy,1);zeros(nZeros_gxy,1)];
mat_gxy = GetAsMatrix(g_vec,n,n);
th_gxy = GetAsVector(GetWithThetas(mat_gxy,n,th1,th2));
th_gxy = th_gxy(1:nCoefficients_gxy);
th_gxy = diag(th_gxy);



% Build the matrix P = DGQ
DPG = D*[ G1*th_fxy*Q1 alpha.*G2*th_gxy*Q2];

end
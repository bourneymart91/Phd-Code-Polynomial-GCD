function Pt = BuildP_BothDegree_STLN(m, m1, m2, n, n1, n2, t, t1, t2, opt_col_index)
% BuildPt(m,m1,m2,n,n1,n2,opt_col,t1,t2)
%
% Build the matrix P_{t}, such that the matrix vector product P*[f;g] gives
% the column c_{t}.
%
% P_{t} * [f;g] = c_{t}
%
% Inputs
%
% m  : (Int) Total degree of polynomial f(x,y)
%
% m1 : (Int) Degree of polynomial f(x,y) with respect to x
%
% m2 : (Int) Degree of polynomial f(x,y) with respect to y
%
% n  : (Int) Total degree of polynomial g(x,y)
%
% n1 : (Int) Degree of polynomial g(x,y) with respect to x
%
% n2 : (Int) Degree of polynomial g(x,y) with respect to y
%
% opt_col : (Int) index of optimal column for removal from Sylvester subresultant
% matrix 
%
% t : (Int) Total degree of polynomial d(x,y)
%
% t1 : (Int) Degree of polynomial d(x,y) with respect to x
%
% t2 : (Int) Degree of polynomial d(x,y) with respect to x



% Get the number of nonZero coefficients in f(x,y), g(x,y) and v(x,y)
nNonZeros_fxy = GetNumNonZeros(m1,m2,m);
nNonZeros_gxy = GetNumNonZeros(n1,n2,n);
nNonZeros_vxy = GetNumNonZeros(n1-t1,n2-t2,n-t);


% Number of columns in T1 of the sylvester matrix
nColumnsT1 = nNonZeros_vxy;

nNonZeros_fv = GetNumNonZeros(m1+n1-t1, m2+n2-t2, m+n-t);

if opt_col_index <= nColumnsT1
    % Optimal column in first partition
    
    % % Build the matrix P
    
    % Build the matrix P1
    P1 = BuildP1_BothDegree_STLN(m, m1, m2, n, n1, n2, t, t1, t2, opt_col_index);
    
    % Build the matrix P2
    P2 = zeros(nNonZeros_fv, nNonZeros_gxy);
    
    
else
    % Optimal column in second partition
    
    % Build the matrix P1
    P1 = zeros(nNonZeros_fv, nNonZeros_fxy);
    
    % Get the position of the optimal column with respect to T(g), so
    % remove the number of columns in T(f).
    opt_col_rel = opt_col_index - nColumnsT1;
    
    % Build the matrix P2
    P2 = BuildP1_BothDegree_STLN(n, n1, n2, m, m1, m2, t, t1, t2, opt_col_rel);
    
    
    
    % Build the matrix P.
    
    
end

Pt = [P1 P2];

end

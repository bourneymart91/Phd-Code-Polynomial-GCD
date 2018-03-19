function Pt = BuildP_RelativeDegree_STLN(m1,m2,n1,n2,idx_opt_col,k1,k2)
% BuildPt(m,m1,m2,n,n1,n2,opt_col,t1,t2)
%
% Build the matrix P_{t}, such that the matrix vector product P*[f;g] gives
% the column c_{t}.
%
% P_{t} * [f;g] = c_{t}
%
% Inputs
%
% m1 : Degree of polynomial f(x,y) with respect to x
%
% m2 : Degree of polynomial f(x,y) with respect to y
%
% n1 : Degree of polynomial g(x,y) with respect to x
%
% n2 : Degree of polynomial g(x,y) with respect to y
%
% idx_opt_col :
%
% k1 : Degree of polynomial d(x,y) with respect to x
%
% k2 : Degree of polynomial d(x,y) with respect to y


% Get the number of coefficients in polynomial f
nCoeff_f = (m1+1).*(m2+1);

% Get the number of coefficients in polynomial g
nCoeff_g = (n1+1).*(n2+1);

% Number of columns in T1 of the sylvester matrix
nColumnsT1 = (n1-k1+1) * (n2-k2+1);

nRows = (m1+n1-k1+1)*(m2+n2-k2+1);

if idx_opt_col <= nColumnsT1
    % Optimal column in first partition
    
    % % Build the matrix P
    
    % Build the matrix P1
    P1 = BuildP1_RelativeDegree_STLN(m1,m2,n1,n2,idx_opt_col,k1,k2);
    
    % Build the matrix P2
    P2 = zeros(nRows,nCoeff_g);
    
    
    
else
    % Optimal column in second partition
    
    % Build the matrix P1
    P1 = zeros(nRows,nCoeff_f);
    
    % Build the matrix P2
    % Get the position of the optimal column with respect to T(g)
    opt_col_rel = idx_opt_col - nColumnsT1;
    P2 = BuildP1_RelativeDegree_STLN(n1,n2,m1,m2,opt_col_rel,k1,k2);
    
    
end

Pt = [P1 P2];

end
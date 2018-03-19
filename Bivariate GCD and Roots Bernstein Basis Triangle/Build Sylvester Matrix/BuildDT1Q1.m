function DT1Q1 = BuildDT1Q1(fxy, m, n_k)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomail f(x,y)
%
% m : (Int) Total degree of f(x,y)
%
% n_k : Degree of v_{k}(x,y), and subscript of T_{n-k}(f(x,y)) 

% Build the diagonal matrix D_{-1}_{m+n-k}
D = BuildD_2Polys(m,n_k);

% Build the matrix T_{n-k}
T1 = BuildT1(fxy,m,n_k);

% Build the matrix Q_{n-k}
Q1 = BuildQ1(n_k);

DT1Q1 = D*T1*Q1;

end
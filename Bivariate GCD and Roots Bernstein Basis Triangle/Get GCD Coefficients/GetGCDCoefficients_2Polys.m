function dxy = GetGCDCoefficients_2Polys(fxy, gxy, uxy, vxy, m, n, t)
% Compute coefficients of polynomial d(x,y), by solving the system Ax = b
% where A is the coefficient matrix given by H^{-1}[C(u);C(v)]G, b is the
% vector of coefficients of f(x,y) and g(x,y) and b is an unknown vector of
% coefficients of d(x,y)
%
% % Inputs.
%
% [fxy, gxy] : Coefficients of polynomial f(x,y) and g(x,y)
% 
% [uxy, vxy] : Coefficients of polynomial u(x,y) and v(x,y)
%
% [m, n] : Total degree of f(x,y) and g(x,y)
%
% t : Total degree of d(x,y)
%
%
% % Outputs.
%
% dxy : Coefficients of polynomail d(x,y)



% % 
% Build the matrix HCG

% Build matrix H^{-1}
H1 = BuildD_2Polys(m-t, t);
H2 = BuildD_2Polys(n-t, t);
H = blkdiag(H1,H2);

% Build matrix C_{t}(u)
C1 = BuildT1(uxy, m-t, t);

% Build matrix C_{t}(v)
C2 = BuildT1(vxy, n-t, t);

C = [C1 ; C2];

% Build matrix G of binomial corresponding to d(x,y)
G = BuildQ1(t);

% %
% Build RHS vector consisting of coefficients of f(x,y) and g(x,y) (Note
% the zeros in the lower right of the matrices of f(x,y) and g(x,y) are 
% first removed)

% Get vector of coefficients of f(x,y)
v_fxy = GetAsVector(fxy);
nCoeffs_fxy = nchoosek(m+2,2);
v_fxy = v_fxy(1:nCoeffs_fxy);

% Get vector of coefficients of g(x,y)
v_gxy = GetAsVector(gxy);
nCoeffs_gxy = nchoosek(n+2,2);
v_gxy = v_gxy(1:nCoeffs_gxy);

% Build RHS Vector
rhs_vec = [v_fxy ; v_gxy];

% % Solve System
% Build System
x_ls = SolveAx_b(H*C*G, rhs_vec);

% % Get coefficients of d(x,y)

% Get number of zeros in d(x,y)
try
    nZeros_dxy = nchoosek(t+1, 2);
catch
    nZeros_dxy = 0;
end

% Append zeros to vector of coefficeints of d(x,y) so that a matrix can be 
% formed
v_dxy = [...
    x_ls;
    zeros(nZeros_dxy, 1)
    ];

% Get matrix d(x,y)
dxy = GetAsMatrix(v_dxy, t, t);

end
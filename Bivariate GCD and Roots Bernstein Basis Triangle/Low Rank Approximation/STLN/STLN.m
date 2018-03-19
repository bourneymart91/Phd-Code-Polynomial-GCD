function [fxy_lr, gxy_lr, uxy_lr, vxy_lr] = STLN(fxy, gxy, k)
% Compute the low rank approximation of the Sylvester matrix S_{k}(f,g)
% by method of Structured Total Least Norm STLN.
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% k : (Int) Degree of common divisor d(x,y) and index of the Sylvester matrix
%     S_{k} whose low rank approximation is considered.
%
% % Outputs.
%
% fxy_lr : (Matrix) Coefficients f(x,y) obtained from the low rank approx of S_{k}(f,g)
%
% gxy_lr : (Matrix) Coefficients g(x,y) obtained from the low rank approx of S_{k}(f,g)
%
% uxy_lr : (Matrix) Coefficients u(x,y) obtained from the low rank approx of S_{k}(f,g)
%
% vxy_lr : (Matrix) Coefficients v(x,y) obtained from the low rank approx of S_{k}(f,g)


% Get the degree of f(x,y) and g(x,y)
[m,~] = GetDegree_Bivariate(fxy);
[n,~] = GetDegree_Bivariate(gxy);

% Get number of coefficients in f(x,y)
nCoefficients_fxy = nchoosek(m + 2,2);

% Get the number of zeros in the matrix containing f(x,y)
nZeros_fxy = nchoosek(m + 1,2);

% Get number of coefficients in g(x,y)
nCoefficients_gxy = nchoosek(n + 2,2);

% Get the number of zeros in the matrix containing g(x,y)
nZeros_gxy = nchoosek(n + 1,2);


% % 
% Build matrices 



% Build the kth Sylvester Subresultant matrix S_{k}(f,g) = D*T(f,g)*Q

% Build the diagonal matrix D^{-1}
D = BuildD_2Polys(m, n - k);

% Build the matrix T_{n-k}(f)
T1_fxy = BuildT1(fxy, m, n - k);

% Build the matrix T_{m-k}(g)
T1_gxy = BuildT1(gxy, n, m - k);

% Build the matrix Q_{k}
Q = BuildQ_2Polys(m,n,k);

% Build the kth Sylvester subresultant matrix.
DTQ_fg = D*[T1_fxy T1_gxy] * Q;

% Get index of optimal column for removal from S_{k}(f,g)
idx_col = GetOptimalColumn(DTQ_fg);

% Remove optimal column from S_{k}(f,g)
Ak_fg = DTQ_fg;
Ak_fg(:, idx_col) = [];

% Get c_{k} optimal column removed from S_{k}(f,g)
ck = DTQ_fg(:, idx_col);


% %

% Get least squares solution x_ls
xk = SolveAx_b(Ak_fg, ck);

% Get the residual vector associated with LS solution
t = (ck) - (Ak_fg * xk);

% %

% Build D*P*Q - The product DPQ * [f;g] returns the column of the Sylvester
% matrix S(f,g).
DPQ = BuildDPQ_STLN(m, n, k, idx_col);

% Get coefficients of polynomial f(x,y)
v_fxy = GetAsVector(fxy);
v_fxy = v_fxy(1 : nCoefficients_fxy);

% Get coefficients of polynomial g(x,y)
v_gxy = GetAsVector(gxy);
v_gxy = v_gxy(1 : nCoefficients_gxy);

test1a = DPQ * [v_fxy; v_gxy];
test1b = ck;
test1 = norm(test1a - test1b);
display(test1);

% Obtain the vector \hat{x} which contains x_ls with a zero inserted into
% the position of the optimal column.
x = [xk(1:idx_col - 1) ; 0 ; xk(idx_col : end)];

% Build matrix Y_{k}
DYQ = BuildDYQ_STLN(x,m,n,k);

test2a = DYQ * [v_fxy; v_gxy];
test2b = DTQ_fg * x;
test2 = norm(test2a - test2b);
display(test2);


% Get matrix of coefficients of polynomial z_{f}(x,y)
zf = zeros(nCoefficients_fxy,1);
vec_zf = [zf ; zeros(nZeros_fxy,1)];
zf_xy = GetAsMatrix(vec_zf,m,m);

% Get matrix of coefficients of polynomial z_{g}(x,y)
zg = zeros(nCoefficients_gxy,1);
vec_zg = [zg ; zeros(nZeros_gxy,1)];
zg_xy = GetAsMatrix(vec_zg,n,n);

% Get vector z = [z_{f} z_{g}]
vec_z = [zf ; zg];

% Build the matrix T_{n-k}(z_{f})
T1_zf = BuildT1(zf_xy, m, n - k);

% Build the matrix T_{m-k}(z_{g})
T1_zg = BuildT1(zg_xy, n, m - k);

% Build the matrix S_{k}(z_{f},z_{g})
Sk_zfzg = D * [T1_zf T1_zg] * Q;
Ak_zfzg = Sk_zfzg;
Ak_zfzg(:, idx_col) = [];

% Build matrices H_{z} and H_{x}
H_z = DYQ - DPQ;
H_x = Ak_fg + Ak_zfzg;

% Build matrix C.
C = [H_z H_x];


% Build the matrix E for LSE Problem
nEntries_x = nchoosek(m - k + 2, 2) + nchoosek(n - k + 2, 2) - 1;
nEntries_z = nchoosek(m + 2, 2) + nchoosek(n + 2, 2);

H1 = eye(nCoefficients_fxy) .* (nchoosek(n - k + 2, 2));
H2 = eye(nCoefficients_gxy) .* (nchoosek(m - k + 2, 2));

E = [blkdiag(H1, H2) zeros(nCoefficients_fxy + nCoefficients_gxy, nEntries_x)];


% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    vec_z;...
    xk;
    ];

% Set yy to be the vector which stores all cummulative perturbations.
yy = start_point;

% Set the initial value of vector f to be zero
s = E * (start_point - yy);

% Initialise the iteration counter
ite = 1;

% Set the termination criterion
vCondition(ite) = norm(t)/norm(ck) ;


global SETTINGS
while vCondition(ite) > SETTINGS.STLN_MAX_ERROR && ite < SETTINGS.STLN_MAX_ITERATIONS
    
    % Increment interation counter
    ite = ite + 1;
    
    % Get small petrubations by LSE
    y_lse = LSE_new(E, s, C, t);
    
    % Increment cummulative peturbations
    yy = yy + y_lse;
    
    % update z and x_ls
    delta_z = y_lse(1 : nEntries_z);
    delta_xk = y_lse(nEntries_z + 1 : end);
    
    % Get z1 and z2
    delta_zf = delta_z(1 : nCoefficients_fxy);
    delta_zg = delta_z(nCoefficients_fxy + 1 : nCoefficients_fxy + nCoefficients_gxy);
    
    % Update vectors z, z_{f} and z_{g}
    vec_z = vec_z + delta_z;
    zf =  zf + delta_zf;
    zg =  zg + delta_zg;
    
    % % Build matrix S_{t}(z1,z2)
    vec_zf = [zf ; zeros(nZeros_fxy, 1)];
    zf_xy = GetAsMatrix(vec_zf, m, m);
    
    vec_zg = [zg ; zeros(nZeros_gxy, 1)];
    zg_xy = GetAsMatrix(vec_zg, n, n);
    
    T1_zf = BuildT1(zf_xy, m, n - k);
    T1_zg = BuildT1(zg_xy, n, m - k);
    
    % Build Sylvester matrix S_{t}(zf,zg)
    Sk_zfzg = D*[T1_zf T1_zg] * Q;
    
    % Get the matrix A_{t}(zf,zg)
    Ak_zfzg = Sk_zfzg;
    Ak_zfzg(:,idx_col) = [];
    
    % Get the vector h_{t}
    hk = Sk_zfzg(:,idx_col);
    
    % Get x_{1} and x_{2}
    
    % Obtain the vector \hat{x} which contains x_ls with a zero inserted into
    % the position of the optimal column.
    xk = xk + delta_xk;
    
    A = Ak_fg + Ak_zfzg;
    b = ck + hk;
    %xk = SolveAx_b(A,b);
    
    %
    x = [xk(1 : idx_col - 1) ; 0 ; xk(idx_col : end)];
    
    % Build the matrix Y_{k}
    DYQ = BuildDYQ_STLN(x, m, n, k);
    
    % % % Build the matrix C in the LSE Problem.
    
    % Build H_z
    H_z = DYQ - DPQ;
    H_x = (Ak_fg + Ak_zfzg);
    
    % Build C
    C = [H_z H_x];
    
    % Get the residual vector
    t = (ck + hk) - ((Ak_fg + Ak_zfzg)*xk);
    
    
    % Update f - used in LSE Problem.
    s = E * (start_point - yy);
    
    % Update the termination criterion
    vCondition(ite) = norm(t)./ norm (ck + hk) ;
    
end


LineBreakLarge()
fprintf([mfilename ' : ' sprintf('STLN Number of iterations required : %i \n',ite)])
LineBreakLarge()
SETTINGS.LOW_RANK_APPROX_REQ_ITE = ite;


% Update z and x_ls
vec_z = yy(1:nEntries_z);

% Get z1 and z2
zf = vec_z(1 : nCoefficients_fxy);
zg = vec_z(nCoefficients_fxy + 1 : end);

% Get Coefficients of polynomials f(x,y) and g(x,y)
vec_zf = [zf ; zeros(nZeros_fxy, 1)];
vec_zg = [zg ; zeros(nZeros_gxy, 1)];

zf_xy = GetAsMatrix(vec_zf, m, m);
zg_xy = GetAsMatrix(vec_zg, n, n);


fxy_lr = fxy + zf_xy;
gxy_lr = gxy + zg_xy;

% Get coefficients of polynomials xu(x,y) and xv(x,y)
% split the vector x into \hat{x}_{1} and \hat{x}_{2}
% Get number of coefficients in x1

% split the vector x into \hat{x}_{1} and \hat{x}_{2}
% Get number of coefficients in x1
nCoefficients_vxy = nchoosek(n - k + 2, 2);
nCoefficients_uxy = nchoosek(m - k + 2, 2);

vec_xvxu = [xk(1 : idx_col - 1) ; -1 ; xk(idx_col : end)];

xv = vec_xvxu(1 : nCoefficients_vxy);
xu = -1.* vec_xvxu(nCoefficients_vxy + 1 : nCoefficients_vxy + nCoefficients_uxy);

% Get x1 as a matrix of coefficients for input into BuildT1() function
nZeros_xv = nchoosek(n - k + 1, 2);
nZeros_xu = nchoosek(m - k + 1, 2);

% Get vectors of coefficients of v and u
vec_xv = [ xv ; zeros(nZeros_xv,1)];
vec_xu = [ xu ; zeros(nZeros_xu,1)];

% Get as matrix
mat_xv = GetAsMatrix(vec_xv, n - k, n - k);
mat_xu = GetAsMatrix(vec_xu, m - k, m - k);

vxy_lr = mat_xv;
uxy_lr = mat_xu;

plotgraphs_STLN()



end




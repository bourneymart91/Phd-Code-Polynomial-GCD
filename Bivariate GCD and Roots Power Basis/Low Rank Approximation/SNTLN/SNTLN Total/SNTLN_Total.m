function [ fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr,th1_lr,th2_lr] = ...
    SNTLN_Total(fxy, gxy, m, n, i_alpha, i_th1, i_th2,k,idx_col)
% Obtain the low rank approximation of the Sylvester matrix D*T_{t}(f,g)*Q =
% S_{t}(f,g)
%
% SNTLN( fxy_matrix,gxy_matrix, i_alpha, i_th1, i_th2,t1,t2,opt_col)
%
% Inputs:
%
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) in the Bernstein basis.
%
% gxy : (Matrix) Coefficients of polynomial g, in the Bernstein basis.
%
% i_alpha : (Float) Initial value of alpha
%
% i_th1 : (Float) Initial value of theta1
%
% i_th2 : (Float) Initial value of theta2
%
% t1 : (Int) Degree of AGCD d(x,y) with respect to x
%
% t2 : (Int) Degree of AGCD d(x,y) with repsect to y
%
% opt_col : (Float) Optimal column for removal from the sylvester matrix, such that col
%           is the column which is most likely a linear combination of the others.
%
%
% % Outputs:
%
%
% fxy_lr : (Matrix) Coefficients of f(x,y) on output, in standard bernstein basis,
% including added structured perturbations.
%
% gxy_lr : (Matrix) Coefficients of g(x,y) on output, in standard bernstein basis,
% including added structured perturbations.
%
% uxy_lr : (Matrix) Coefficients of the polynomial u(x,y)
%
% vxy_lr : (Matrix) Coefficients of the polynomial v(x,y)
%
% alpha_lr : (Float) Optimal value of \alpha
%
% th1_lr : (Float) Optimal value of \theta_{1}
%
% th2_lr : (Float) Optimal value of \theta_{2}


% Global Inputs
global SETTINGS

% %
% Pad the coefficients of fxy and gxy
% this is equivalent to degree elevating so that f is of degree (m,m), and
% g is of degree (n,n)
fxy_matrix_padd = zeros(m+1,m+1);
gxy_matrix_padd = zeros(n+1,n+1);


[nRows_fxy,nCols_fxy] = size(fxy);
fxy_matrix_padd(1:nRows_fxy,1:nCols_fxy) = fxy;

[nRows_gxy,nCols_gxy] = size(gxy);
gxy_matrix_padd(1:nRows_gxy,1:nCols_gxy) = gxy;

fxy = fxy_matrix_padd;
gxy = gxy_matrix_padd;

% %

% Set the initial iterations number
ite = 1;

% Set initial values of alpha and theta
th1(ite) = i_th1;
th2(ite) = i_th2;
alpha(ite) = i_alpha;

% Get the number of coefficients in the polynomial f(x,y)
nNonZeros_fxy = nchoosek(m+2,2);
nZeros_fxy = nchoosek(m+1,2);

% Get the number of coefficients in the polynomial g(x,y)
nNonZeros_gxy = nchoosek(n+2,2);
nZeros_gxy = nchoosek(n+1,2);

% Get the number of coefficients in both f(x,y) and g(x,y)
nNonZeros_fg = nNonZeros_fxy + nNonZeros_gxy;

% Get the number of coefficients in v(x,y)
nNonZeros_vxy = nchoosek(n-k+2,2);
nZeros_vxy = nchoosek(n-k+1,2);

% Get the number of coefficients in u(x,y)
nNonZeros_uxy = nchoosek(m-k+2,2);
nZeros_uxy = nchoosek(m-k+1,2);

% Get the number of coefficients in the unknown vector x, where A_{t}x =
% c_{t}.
nNonZeros_x = nNonZeros_uxy + nNonZeros_vxy - 1;

% Create the identity matrix I, the matrix M formed from I by removing the
% column equivalent to the optimal column for removal from the Sylvester
% subresultant matrix, so that S_{t}(f,g)*M = A_{t}, where A_{t} is the
% Sylvester subresultant matrix with the column removed.

% Get the number of columns in C_{t}(f), the first partition of the Sylvester
% Matrix S_{t}(f,g)
nColumns_T1 = nNonZeros_vxy;

% Get the number of columns in C_{t}(g), the second partition of the
% Sylvester matrix S_{t}(f,g)
nColumns_T2 = nNonZeros_uxy;

% Get the total number of columns in the Sylvester matrix S_{k}(f,g)
nColumns_Sk = nColumns_T1 + nColumns_T2;

% Get the number of rows in the Sylvester subresultant matrix S_{k}(f,g)
nRows_Sk = nchoosek(m+n-k+2,2);


% Create the identity matrix
I = eye(nColumns_Sk, nColumns_Sk);

% Create the matrix M, such that S(f,g)*M gives A_{t}, the Sylvester Matrix
% with the optimal column removed.
M = I;
M(:,idx_col) = [];

% Let e be the column removed from the identity matrix, such that
% S_{t}(f,g) * e gives the column c_{t}, where c_{t} is the optimal column
% removed from the Sylvester subresultant.
e = I(:,idx_col);

% % Preprocessing
% Obtain polynomials in Modified Bernstein Basis, using initial values of
% alpha and theta.

% Multiply the rows of fxy_matrix by theta1, and multiply the cols of
% fxy_matrix by theta2
fww_matrix = GetWithThetas(fxy, th1(ite), th2(ite));
gww_matrix = GetWithThetas(gxy, th1(ite), th2(ite));

% Form the Coefficient Matrix T = [C(f(w1,w2))|alpha * C(g(w1,w2))] such that T*x = [col]
T1_fxy = BuildT1_Total_Bivariate(fww_matrix, m, n-k);
T2_gxy = BuildT1_Total_Bivariate(gww_matrix, n, m-k);
T_fg = [T1_fxy alpha(ite).*T2_gxy];

%
% Calculate the partial derivatives of f(w,w) and g(w,w) with respect to \alpha
Partial_fw_wrt_alpha            = zeros(m+1,m+1);
Partial_alpha_gw_wrt_alpha      = gxy;

%
% Calculate the partial derivatives of f(w,w) with respect to \theta_1
theta_mat = diag((0:1:m) ./ th1(ite));
Partial_fw_wrt_th1    = theta_mat * fww_matrix;

% Calculate the partial derivative of g(w,w) with respect to theta_1
theta_mat = diag((0:1:n) ./ th1(ite));
Partial_gw_wrt_th1    = theta_mat * gww_matrix;

% Calculate the partial derivative of f(w,w) with respect to theta_2
theta_mat = diag((0:1:m) ./ th2(ite));
Partial_fw_wrt_th2 = fww_matrix * theta_mat;


% Calculate the partial deriviates of g(w,w) with respect to theta_2
theta_mat = diag((0:1:n)./ th2(ite));
Partial_gw_wrt_th2 = gww_matrix * theta_mat;

% Build the derivative of T(f,g) with respect to alpha
T1_f_wrt_alpha = BuildT1_Total_Bivariate(Partial_fw_wrt_alpha, m, n-k);
T2_g_wrt_alpha = BuildT1_Total_Bivariate(Partial_alpha_gw_wrt_alpha, n, m-k);
T_fg_wrt_alpha = [T1_f_wrt_alpha T2_g_wrt_alpha];

% Calculate the derivative of T(f,g) with respect to theta_{1}
T1_f_wrt_th1 = BuildT1_Total_Bivariate(Partial_fw_wrt_th1, m, n-k);
T2_g_wrt_th1 = BuildT1_Total_Bivariate(Partial_gw_wrt_th1, n, m-k);
Partial_T_fg_wrt_th1 = [T1_f_wrt_th1 alpha(ite)* T2_g_wrt_th1];

% Calcualte the derivative of T(f,g) with respect to theta_2
T1_f_wrt_th2 = BuildT1_Total_Bivariate(Partial_fw_wrt_th2, m, n-k);
T2_g_wrt_th2 = BuildT1_Total_Bivariate(Partial_gw_wrt_th2, n, m-k);
Partial_T_fg_wrt_th2 = [T1_f_wrt_th2 alpha(ite)*T2_g_wrt_th2];

% Initialise the vector z of structured perturbations
% if we are working with strictly the roots problem, the number of entries
% in z can be reduced.
zk = zeros(nNonZeros_fg , 1);


% Initilaise the derivative of N wrt alpha.
Partial_N_wrt_alpha = zeros(nRows_Sk, nColumns_T1 + nColumns_T2);

% Initilaise the derivative of N wrt theta_1.
Partial_N_wrt_th1 = zeros(nRows_Sk, nColumns_T1 + nColumns_T2);

% Initialise the derivative of N wrt theta 2
Partial_N_wrt_th2 = zeros(nRows_Sk, nColumns_T1 + nColumns_T2);

%
% Initialise the derivative of h
% Calculate the derivatives wrt alpha and theta of the column of DNQ
% that is moved to the right hand side.
Partial_h_wrt_alpha = Partial_N_wrt_alpha*e;
Partial_h_wrt_theta1 = Partial_N_wrt_th1*e;
Partial_h_wrt_theta2 = Partial_N_wrt_th2*e;


% Get the matrix A_{k}(f,g), which is the subresultant matrix S(f,g) with
% an opitmal column removed
Ak = T_fg;
ck = T_fg(:,idx_col);
Ak(:,idx_col) = [];

% Build the matrix P
Pk = BuildP_TotalDegree_SNTLN(m, n, k, alpha(ite), th1(ite), th2(ite), idx_col);

% Test 1 - Test the Build P function
% Get the coefficients of f(x,y) in matrix form
fxy_vec = GetAsVector_Version1(fxy);
fxy_vec = fxy_vec(1:nNonZeros_fxy);
gxy_vec = GetAsVector_Version1(gxy);
gxy_vec = gxy_vec(1:nNonZeros_gxy);

%Test
test1a = ck;
test1b = Pk*[fxy_vec;gxy_vec];
test1 = norm(test1a - test1b);
display(test1)


% Calculate the derivatives of ck wrt alpha and theta.
Partial_ck_wrt_alpha = T_fg_wrt_alpha*e;
Partial_ck_wrt_th1 = Partial_T_fg_wrt_th1*e;
Partial_ck_wrt_th2 = Partial_T_fg_wrt_th2*e;

% Perform QR decomposition of Ak to obtain the solution x
xk = SolveAx_b(Ak, ck);

% Get the vector x such that S_{k}(f,g) * x = ck
first_part = xk(1:(idx_col-1));
second_part = xk(idx_col:end);
x = [first_part ; 0 ; second_part];

% Build Matrix Y, where Y(x1,x2)*[f;g] = S(f,g)*[x1;x2]
Yk = BuildY_TotalDegree_SNTLN(x, m, n, k, alpha(ite), th1(ite), th2(ite));

% Test Y
test2a  = Yk * [fxy_vec;gxy_vec];
test2b = T_fg * x;
test2 = norm(test2a - test2b);
display(test2)

% Calculate the initial residual r = ck - (Ak*x)
res_vec = ck - (T_fg*M*xk);

% Get the matrix p, which will store all the perturbations returned from LSE file
nEntries = nNonZeros_fxy...
    + nNonZeros_gxy ...
    + nNonZeros_x ...
    + 3;


% Set the initial value of vector p to be zero
f = zeros(nNonZeros_fxy...
    + nNonZeros_gxy ...
    + nNonZeros_x ...
    + 3,1);

%
% Set the intial value of E to the identity matrix
E = eye(nEntries);

%
% Create the matrix (T+N), initially N is empty so this is the same as T.
T1_fxy = BuildT1_Total_Bivariate(fww_matrix, m, n-k);
T2_gxy = BuildT1_Total_Bivariate(gww_matrix, n, m-k);
TN = [T1_fxy alpha(ite).*T2_gxy];

%
% Create The matrix (T+N) with respect to alpha
T1_f_wrt_alpha = BuildT1_Total_Bivariate(Partial_fw_wrt_alpha, m, n-k);
T2_g_wrt_alpha = BuildT1_Total_Bivariate(Partial_alpha_gw_wrt_alpha, n, m-k);
TN_wrt_alpha = [T1_f_wrt_alpha T2_g_wrt_alpha];

%
% Create The matrix (T+N) with respect to theta1
T1_f_wrt_th1 = BuildT1_Total_Bivariate(Partial_fw_wrt_th1, m, n-k);
T2_g_wrt_th1 = BuildT1_Total_Bivariate(Partial_gw_wrt_th1, n, m-k);
TN_wrt_theta1 = [T1_f_wrt_th1 alpha(ite) * T2_g_wrt_th1];

%
% Create The matrix (T+N) with respect to theta2
T1_f_wrt_th2 = BuildT1_Total_Bivariate(Partial_fw_wrt_th2, m, n-k);
T2_g_wrt_th2 = BuildT1_Total_Bivariate(Partial_gw_wrt_th2, n, m-k);
TN_wrt_theta2 = [T1_f_wrt_th2 alpha(ite) * T2_g_wrt_th2];

%
% Create the matrix C for input into iteration

H_z     = Yk-Pk;

H_x     = TN*M;

H_alpha  = TN_wrt_alpha*M*xk - ...
    (Partial_ck_wrt_alpha + Partial_h_wrt_alpha);

H_theta1 = TN_wrt_theta1*M*xk - ...
    (Partial_ck_wrt_th1 + Partial_h_wrt_theta1);

H_theta2 = TN_wrt_theta2*M*xk - ...
    (Partial_ck_wrt_th2 + Partial_h_wrt_theta2);

C       = [H_z H_x H_alpha H_theta1 H_theta2];

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    zk;...
    xk;...
    alpha(ite);...
    th1(ite);...
    th2(ite)
    ];

yy =  start_point;

% Set the termination criterion to a large value. It will be
% over written later.
condition(ite) = norm(res_vec) ./ norm(ck);



while condition(ite) >(SETTINGS.MAX_ERROR_SNTLN) &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN
       
    
    % Use the QR decomposition to solve the LSE problem
    % min |y-p| subject to Cy=q
    y = LSE(E, f, C, res_vec);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % Break down y into its sections
    
    % Get the coefficients corresponding to f and g
    delta_zk = y(1:nNonZeros_fxy + nNonZeros_gxy ,1);
    
    % Remove the zk coefficients from the list of coefficients
    y(1:nNonZeros_fxy + nNonZeros_gxy) = [];
    
    % Get the coefficients corresponding to x
    delta_xk = y(1:nNonZeros_x,1);
    
    % Remove them from the list of coefficients
    y(1:nNonZeros_x) = [];
    
    % Get the coefficient corresponding to \alpha
    delta_alpha = y(1:1);
    y(1) = [];
    
    % Get the coefficient corresponding to \theta_{1}
    delta_th1    = y(1:1);
    y(1) = [];
    
    % Get the coefficient corresponding to \theta_{2}
    delta_th2    = y(1:1);
    y(1) = [];
    
    % % Update the variables
    
    % Update variables z_{k}, where z_{k} are perturbations in the
    % coefficients of f and g.
    zk = zk + delta_zk;
    
    % Update x_{k}, where x_{k} is the solution vector, containing
    % coefficients u and v.
    xk = xk + delta_xk;
    
    % Update alpha
    alpha(ite) = alpha(ite-1) + delta_alpha;
    
    % Update theta_{1}
    th1(ite) = th1(ite-1) + delta_th1;
    
    % Update theta_{2}
    th2(ite) = th2(ite-1) + delta_th2;
    
    % Obtain polynomials in modified bersntein basis a_{i}\theta^{i}
    
    % Obtain new f(w,w) with improved theta1, and theta2
    fww_matrix = GetWithThetas(fxy, th1(ite), th2(ite));
    
    % Obtain new g(w,w) with improved theta1 and theta2
    gww_matrix = GetWithThetas(gxy, th1(ite), th2(ite));
    
    % Construct the Sylvester subresultant matrix S.
    T1_fxy = BuildT1_Total_Bivariate(fww_matrix, m, n-k);
    T2_gxy = BuildT1_Total_Bivariate(gww_matrix, n, m-k);
    T_fg = [T1_fxy alpha(ite).*T2_gxy];
    
    % Calculate the partial derivatives of fw and gw with respect to alpha
    Partial_fw_wrt_alpha            = zeros(m+1,m+1);
    Partial_alpha_gw_wrt_alpha      = gww_matrix;
    
    % %
    % Calculate the partial derivatives of fw and gw with respect to theta1
    % divide the rows by theta1 and multiply by the old power
    
    % Get the partial derivative of f with respect to \theta_{1}
    temp_mat = diag((0:1:m)./th1(ite));
    Partial_fw_wrt_th1 = temp_mat * fww_matrix;
    
    % Get the partial derivative of g with respect to \theta_{1}
    temp_mat = diag((0:1:n)./th1(ite));
    Partial_gw_wrt_th1 = temp_mat * gww_matrix;
    
    % Get the partial derivative of f with respect to \theta_{2}
    temp_mat = diag((0:1:m)./th2(ite));
    Partial_fw_wrt_th2 =  fww_matrix * temp_mat;
    
    % Get the partial derivative of g with respect to \theta_{2}
    temp_mat = diag((0:1:n)./th2(ite));
    Partial_gw_wrt_th2 =  gww_matrix * temp_mat;
    
    % Calculate the Partial derivative of T with respect to alpha.
    T1_f_wrt_alpha = BuildT1_Total_Bivariate(Partial_fw_wrt_alpha, m, n-k);
    T2_g_wrt_alpha = BuildT1_Total_Bivariate(Partial_alpha_gw_wrt_alpha, n, m-k);
    Partial_T_wrt_alpha = [T1_f_wrt_alpha T2_g_wrt_alpha];
    
    % Calculate the partial derivative of T with respect to theta1
    T1_f_wrt_th1 = BuildT1_Total_Bivariate(Partial_fw_wrt_th1, m, n-k);
    T2_g_wrt_th1 = BuildT1_Total_Bivariate(Partial_gw_wrt_th1, n, m-k);
    Partial_T_fg_wrt_th1 = [T1_f_wrt_th1 alpha(ite)*T2_g_wrt_th1];
    
    % Calculate the partial derivative of T with respect to theta2
    T1_f_wrt_th2 = BuildT1_Total_Bivariate(Partial_fw_wrt_th2, m, n-k);
    T2_g_wrt_th2 = BuildT1_Total_Bivariate(Partial_gw_wrt_th2, n, m-k);
    Partial_T_fg_wrt_th2 = [T1_f_wrt_th2 alpha(ite)*T2_g_wrt_th2];
    
    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ck = T_fg*e;
    
    % Calculate the derivatives of c_{k} with respect to \alpha and \theta
    Partial_ck_wrt_alpha = Partial_T_wrt_alpha*e;
    Partial_ck_wrt_th1 = Partial_T_fg_wrt_th1*e;
    Partial_ck_wrt_th2 = Partial_T_fg_wrt_th2*e;
    
    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    z_fx      = zk(1:nNonZeros_fxy);
    z_gx      = zk(nNonZeros_fxy + 1 :end);
    
    % Get the matrices of z_f(x,y) and z_g(x,y)
    z_fx_mat = GetAsMatrix_Version1([z_fx;zeros(nZeros_fxy,1)], m, m);
    z_gx_mat = GetAsMatrix_Version1([z_gx;zeros(nZeros_gxy,1)], n, n);
    
    % Get the matrices of z_f(w,w) and z_g(w,w)
    z_fw_mat = GetWithThetas(z_fx_mat, th1(ite), th2(ite));
    z_gw_mat = GetWithThetas(z_gx_mat, th1(ite), th2(ite));
    
    % Calculate the derivatives of z_f(w,w) and z_g(w,w) with repect to alpha.
    Partial_zfw_wrt_alpha    = zeros(m+1,m+1);
    Partial_zgw_wrt_alpha    = z_gw_mat;
    
    % Calculate the derivative of z_f(w,w) with respect to \theta_{1}.
    temp_mat = diag((0:1:m)./th1(ite));
    Partial_zfw_wrt_theta1 = temp_mat * z_fw_mat;
    
    % Calculate the derivative of z_f(w,w) with respect to \theta_{2}
    temp_mat = diag((0:1:m)./th2(ite));
    Partial_zfw_wrt_theta2 = z_fw_mat * temp_mat;
    
    % Calculate the derivative of z_g(w,w) with respect ot \theta_{1}
    temp_mat = diag((0:1:n)./th2(ite));
    Partial_zgw_wrt_theta1 = temp_mat * z_gw_mat;
    
    % Calculate the deriviate of z_g(w,w) with respect to \theta_{2}
    temp_mat = diag((0:1:n)./th2(ite));
    Partial_zgw_wrt_theta2 = z_gw_mat * temp_mat;
    
    % Build the coefficient Matrix N = [T(z1) T(z2)], of structured 
    % perturbations, with the same structure as T_fg
    N1 = BuildT1_Total_Bivariate(z_fw_mat, m, n-k);
    N2 = BuildT1_Total_Bivariate(z_gw_mat, n, m-k);
    N = [N1 alpha(ite).*N2];
    
    % Build the coefficient matrix N with respect to alpha
    T1_zf_wrt_alpha = BuildT1_Total_Bivariate(Partial_zfw_wrt_alpha, m, n-k);
    T2_zg_wrt_alpha = BuildT1_Total_Bivariate(Partial_zgw_wrt_alpha, n, m-k);
    Partial_N_wrt_alpha = [T1_zf_wrt_alpha T2_zg_wrt_alpha];
    
    % Calculate the derivatives of DNQ with respect to theta
    T1_zf_wrt_th1 = BuildT1_Total_Bivariate(Partial_zfw_wrt_theta1, m, n-k);
    T2_zg_wrt_th1 = BuildT1_Total_Bivariate(Partial_zgw_wrt_theta1, n, m-k);
    Partial_N_wrt_th1 = [T1_zf_wrt_th1 alpha(ite).*T2_zg_wrt_th1];
    
    % Calculate the derivatives of DNQ with respect to theta
    T1_zf_wrt_th2 = BuildT1_Total_Bivariate(Partial_zfw_wrt_theta2, m, n-k);
    T2_zg_wrt_th2 = BuildT1_Total_Bivariate(Partial_zgw_wrt_theta2, n, m-k);
    Partial_N_wrt_th2 = [T1_zf_wrt_th2 alpha(ite).*T2_zg_wrt_th2];
    
    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{k} the column of S_{k} moved to the RHS
    hk = N*e;
    
    % Calculate the derivative of h with respect to alpha
    h_alpha = Partial_N_wrt_alpha*e;
    % Calculate the derivative of h with respect to theta1
    h_theta1 = Partial_N_wrt_th1*e;
    % Calculate the derivative of h with respect to theta2
    h_theta2 = Partial_N_wrt_th2*e;
    
    % Build the matrix (T+N)
    T1_fxy = BuildT1_Total_Bivariate(fww_matrix + z_fw_mat, m, n-k);
    T2_gxy = BuildT1_Total_Bivariate(gww_matrix + z_gw_mat, n, m-k);
    TN = [T1_fxy alpha(ite).*T2_gxy];
    
    % Calculate the paritial derivative of (T+N) with respect to
    % alpha
    TN1_wrt_alpha = BuildT1_Total_Bivariate(Partial_fw_wrt_alpha + Partial_zfw_wrt_alpha, m, n-k);
    TN2_wrt_alpha = BuildT1_Total_Bivariate(Partial_alpha_gw_wrt_alpha + Partial_zgw_wrt_alpha, n, m-k);
    TN_alpha = [TN1_wrt_alpha TN2_wrt_alpha];
    
    
    % Calculate the paritial derivative of (T+N) with respect to
    % theta1
    TN1_wrt_th1 = BuildT1_Total_Bivariate(Partial_fw_wrt_th1 + Partial_zfw_wrt_theta1, m, n-k);
    TN2_wrt_th1 = BuildT1_Total_Bivariate(Partial_gw_wrt_th1 + Partial_zgw_wrt_theta1, n, m-k);
    TN_theta1 = [TN1_wrt_th1 alpha(ite).*TN2_wrt_th1];
    
    % Calculate the paritial derivative of (T+N) with respect to
    % theta2
    TN1_wrt_th2 = BuildT1_Total_Bivariate(Partial_fw_wrt_th2 + Partial_zfw_wrt_theta2, m, n-k);
    TN2_wrt_th2 = BuildT1_Total_Bivariate(Partial_gw_wrt_th2 + Partial_zgw_wrt_theta2, n, m-k);
    TN_theta2 = [TN1_wrt_th2 alpha(ite).*TN2_wrt_th2];
    
    
    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    first_part = xk(1:(idx_col-1));
    second_part = xk(idx_col:end);
    x = [first_part ; 0 ; second_part];
    
    Yk = BuildY_TotalDegree_SNTLN(x, m, n, k, alpha(ite), th1(ite), th2(ite));
    
    % Calculate the matrix P where ck = P * [f,g]
    Pk = BuildP_TotalDegree_SNTLN(m, n, k, alpha(ite), th1(ite), th2(ite), idx_col);
    
    
    % %
    
    % Get residual as a vector
    rk = (ck+hk) - TN*M*xk ;
    
    % Create the matrix C. This is made up of five submatrices, HZ, Hx,
    % H_alpha and H_theta1 and H_theta2.
    
    Hz          = Yk-Pk;
    
    Hx          = TN*M;
    
    H_alpha     = TN_alpha*M*xk - (Partial_ck_wrt_alpha + h_alpha);
    
    H_theta1    = TN_theta1*M*xk - (Partial_ck_wrt_th1 + h_theta1);
    
    H_theta2    = TN_theta2*M*xk - (Partial_ck_wrt_th2 + h_theta2);
    
    C = [Hz,Hx,H_alpha,H_theta1, H_theta2];  % the matrix C
    
    % Update Residual vector which is used in LSE Problem.
    res_vec = rk;
    
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(rk) / norm(ck + hk);
    
    % Update f which is used in LSE Problem.
    f = -(yy-start_point);
    
end

PlotGraphs_SNTLN()
SETTINGS.LOW_RANK_APPROX_REQ_ITE = ite;

%

% Once iterations are complete, assign fx output, gx output, solution X
% output, alpha output and theta output.

% get the vector zk
zPert_f_vec = zk(1:nNonZeros_fxy);
zPert_f_mat = GetAsMatrix_Version1([zPert_f_vec; zeros(nZeros_fxy,1)],m,m);

zPert_g_vec = zk(nNonZeros_fxy+1:end);
zPert_g_mat = GetAsMatrix_Version1([zPert_g_vec; zeros(nZeros_gxy,1)],n,n);

% Set outputs of low rank approximation
fxy_lr = fxy + zPert_f_mat;
gxy_lr = gxy + zPert_g_mat;

% Get coefficients of the polynomials u(x,y) and v(x,y)

first_part = xk(1:(idx_col-1));
second_part = xk(idx_col:end);
x = [first_part ; -1 ; second_part];

%
vec_vww = x(1:nNonZeros_vxy);
vec_uww = -1.* x(nNonZeros_vxy+1:end);

%
vww = GetAsMatrix_Version1([vec_vww ; zeros(nZeros_vxy,1)],n-k,n-k);
vxy_lr = GetWithoutThetas(vww,th1(ite),th2(ite));

%
uww = GetAsMatrix_Version1([vec_uww ; zeros(nZeros_uxy,1)],m-k,m-k); 
uxy_lr = GetWithoutThetas(uww,th1(ite),th2(ite));

% %
%
alpha_lr = alpha(ite);
th1_lr = th1(ite);
th2_lr = th2(ite);


% Print the number of iterations
LineBreakLarge();
fprintf([mfilename ' : ' sprintf('SNTLN converged within %i iterations \n', ite)]);
LineBreakLarge()

end


















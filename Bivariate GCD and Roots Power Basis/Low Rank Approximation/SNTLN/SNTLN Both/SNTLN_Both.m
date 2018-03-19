function [ fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr,th1_lr,th2_lr] = ...
    SNTLN_Both( fxy, gxy, i_alpha, i_th1, i_th2,m,n,k,k1,k2,idx_col)
% Obtain the low rank approximation of the Sylvester matrix D*T_{t}(f,g)*Q =
% S_{t}(f,g)
%
% SNTLN( fxy_matrix,gxy_matrix, i_alpha, i_th1, i_th2,t1,t2,opt_col)
%
% % Inputs:
%
% fxy : (Matrix) Coefficients of polynomial f, in standard bernstein basis.
%
% gxy : (Matrix) Coefficients of polynomial g, in standard bernstein basis.
%
% i_alpha : (Float) Initial value of alpha
%
% i_th1 : (Float) Initial value of theta1
%
% i_th2 : (Float) Initial value of theta2
%
% m : (Int) Total degree of polynomial f(x,y)
%
% n : (Int) Total degree of polynomial g(x,y)
%
% k : (Int) Total degree of polynomial d(x,y)
%
% k1 : (Int) Degree of AGCD d(x,y) with respect to x
%
% k2 : (Int) Degree of AGCD d(x,y) with repsect to y
%
% idx_col : (Int) Optimal column for removal from the sylvester matrix, such that col
%           is the column which is most likely a linear combination of the others.
%
%
% % Outputs:
%
% fxy_lr : (Matrix) Coefficients of f(x,y) on output, in standard bernstein basis,
% including added structured perturbations.
%
% gxy_lr : (Matrix) Coefficients of g(x,y) on output, in standard bernstein basis,
% including added structured perturbations.
%
% uxy_lr : (Matrix) 
%
% vxy_lr : (Matrix)
% 
% o_alpha : (Float) Optimal value of \alpha
%
% o_th1 : (Float) Optimal value of \theta_{1}
%
% o_th2 : (Float) Optimal value of \theta_{2}


% Global Inputs

global SETTINGS

% Set the initial iterations number
ite = 1;

% Set initial values of \alpha, \theta_{1} and \theta_{2}
th1(ite) = i_th1;
th2(ite) = i_th2;
alpha(ite) = i_alpha;

% Get degree of polynomials f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degree of polynomial g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Get the number of coefficients in the polynomial f(x,y)
nCoefficients_fxy =  (m1+1) * (m2+1);
nNonZeros_fxy = GetNumNonZeros(m1,m2,m);
nZeros_fxy = nCoefficients_fxy - nNonZeros_fxy;

% Get the number of coefficients in the polynomial g(x,y)
nCoefficients_gxy = (n1+1) * (n2+1);
nNonZeros_gxy = GetNumNonZeros(n1,n2,n);
nZeros_gxy = nCoefficients_gxy - nNonZeros_gxy;

% Get the number of coefficients of the product f*v or g*u
nCoefficients_fv = GetNumNonZeros(m1+n1-k1,m2+n2-k2,m+n-k);

% Get the number of coefficients in v(x,y)
nCoefficients_vxy = (n1-k1+1) * (n2-k2+1);
nNonZeros_vxy = GetNumNonZeros(n1-k1,n2-k2,n-k);
nZeros_vxy = nCoefficients_vxy - nNonZeros_vxy;

% Get the number of coefficients in u(x,y)
nCoefficients_uxy = (m1-k1+1) * (m2-k2+1);
nNonZeros_uxy = GetNumNonZeros(m1-k1,m2-k2,m-k);
nZeros_uxy = nCoefficients_uxy - nNonZeros_uxy;

% Get the number of coefficients in the unknown vector x, where A_{t}x =
% c_{t}.
nCoefficients_x = nNonZeros_uxy + nNonZeros_vxy - 1;

% Get the number of columns in T_{k}(f), the first partition of the Sylvester
% Matrix S_{k}(f,g)
nColumns_Tf = nNonZeros_vxy;

% Get the number of columns in T_{k}(g), the second partition of the
% Sylvester matrix S_{k}(f,g)
nColumns_Tg = nNonZeros_uxy;

% Get the total number of columns in the Sylvester matrix S_{k}(f,g)
nColumns_Sylv = nColumns_Tf + nColumns_Tg;

nRows_Sylv_mat = nCoefficients_fv;

% Create the identity matrix
I = eye(nColumns_Sylv, nColumns_Sylv);

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

% Multiply the rows of fxy and gxy by theta1, and multiply the cols of
% fxy and gxy by theta2
fww = GetWithThetas(fxy, th1(ite), th2(ite));
gww = GetWithThetas(gxy, th1(ite), th2(ite));

% Form the Coefficient Matrix T = [C(f(w1,w2))|alpha * C(g(w1,w2))] such that T*x = [col]
Tf = BuildT1_Both_Bivariate(fww, m, n-k, n1-k1, n2-k2);
Tg = BuildT1_Both_Bivariate(gww, n, m-k, m1-k1, m2-k2);
Sk_fg = [Tf alpha(ite).*Tg];

%
% Calculate the partial derivatives of f(w,w) and g(w,w) with respect to \alpha
Partial_fww_wrt_alpha            = zeros(m1+1, m2+1);
Partial_alpha_gww_wrt_alpha      = gxy;

%
% Calculate the partial derivatives of f(w,w) with respect to \theta_1
theta_mat = diag((0:1:m1) ./ th1(ite));
Partial_fww_wrt_th1    = theta_mat * fww;

%
% Calculate the partial derivative of g(w,w) with respect to theta_1
theta_mat = diag((0:1:n1) ./ th1(ite));
Partial_gww_wrt_th1    = theta_mat * gww;

%
% Calculate the partial derivative of f(w,w) with respect to theta_2
theta_mat = diag((0:1:m2) ./ th2(ite));
Partial_fww_wrt_th2 = fww * theta_mat;

%
% Calculate the partial deriviates of g(w,w) with respect to theta_2
theta_mat = diag((0:1:n2)./ th2(ite));
Partial_gww_wrt_th2 = gww * theta_mat;

% %
% Build the derivative of T(f,g) with respect to alpha
Tfww_wrt_alpha = BuildT1_Both_Bivariate(Partial_fww_wrt_alpha, m, n-k, n1-k1, n2-k2);
Tgww_wrt_alpha = BuildT1_Both_Bivariate(Partial_alpha_gww_wrt_alpha, n, m-k, m1-k1, m2-k2);
Tfg_alpha = [Tfww_wrt_alpha Tgww_wrt_alpha];

% %
% Calculate the derivative of T(f,g) with respect to theta_{1}
Tfww_wrt_th1 = BuildT1_Both_Bivariate(Partial_fww_wrt_th1, m, n-k, n1-k1 ,n2-k2);
Tgww_wrt_th1 = BuildT1_Both_Bivariate(Partial_gww_wrt_th1, n, m-k, m1-k1, m2-k2);
Partial_Tfg_wrt_th1 = [Tfww_wrt_th1 alpha(ite)* Tgww_wrt_th1];

% %
% Calcualte the derivative of T(f,g) with respect to theta_2
Tfww_wrt_th2 = BuildT1_Both_Bivariate(Partial_fww_wrt_th2, m, n-k, n1-k1, n2-k2);
Tgww_wrt_th2 = BuildT1_Both_Bivariate(Partial_gww_wrt_th2, n, m-k, m1-k1, m2-k2);
Partial_Tfg_wrt_th2 = [Tfww_wrt_th2 alpha(ite)*Tgww_wrt_th2];

% Initialise the vector z of structured perturbations
zk = zeros(nNonZeros_fxy + nNonZeros_gxy , 1);

%
% Initilaise the derivative of N wrt alpha.
Partial_N_wrt_alpha   = zeros(nRows_Sylv_mat,nColumns_Tf + nColumns_Tg);

% Initilaise the derivative of N wrt theta_1.
Partial_N_wrt_th1   = zeros(nRows_Sylv_mat,nColumns_Tf + nColumns_Tg);

% Initialise the derivative of N wrt theta 2
Partial_N_wrt_th2   = zeros(nRows_Sylv_mat,nColumns_Tf + nColumns_Tg);

%
% Initialise the derivative of h
% Calculate the derivatives wrt alpha and theta of the column of DNQ
% that is moved to the right hand side.
Partial_h_wrt_alpha     = Partial_N_wrt_alpha*e;
Partial_h_wrt_theta1    = Partial_N_wrt_th1*e;
Partial_h_wrt_theta2    = Partial_N_wrt_th2*e;


% Get the matrix A_{k}(f,g), which is the subresultant matrix S(f,g) with
% an opitmal column removed
Ak_fg = Sk_fg;
ck = Sk_fg(:,idx_col);
Ak_fg(:,idx_col) = [];

% Build the matrix P
Pk = BuildP_BothDegree_SNTLN(m, m1, m2, n, n1, n2, k, k1, k2, alpha(ite), th1(ite), th2(ite), idx_col);

% Get the coefficients of f(x,y) in matrix form
fxy_vec = GetAsVector_Version1(fxy);
fxy_vec = fxy_vec(1:nNonZeros_fxy);

gxy_vec = GetAsVector_Version1(gxy);
gxy_vec = gxy_vec(1:nNonZeros_gxy);

% Test P
test1_a = ck;
test1_b =  Pk * [fxy_vec; gxy_vec];
test1 = norm(test1_a - test1_b);
display(test1);

% %
% Calculate the derivatives of ck wrt alpha and theta.
Partial_ck_wrt_alpha = Tfg_alpha*e;
Partial_ck_wrt_th1 = Partial_Tfg_wrt_th1*e;
Partial_ck_wrt_th2 = Partial_Tfg_wrt_th2*e;


% %
% Perform QR decomposition of Ak to obtain the solution x
xk = SolveAx_b(Ak_fg,ck);

% Insert a zero into the least squares solution x_ls so that 
% S_{k,k1,k2} x = c_{k,k1,k2}. Also so that when x is split into two 
% vectors x1 and x2. S(x1,x2) [f;g] =  ck.
first_part = xk(1:(idx_col-1));
second_part = xk(idx_col:end);
x = [first_part ; 0 ; second_part];

% Build Matrix Y, where Y(x1,x2)*[f;g] = S(f,g)*[x1;x2]
Yk = BuildY_BothDegree_SNTLN(x, m, m1, m2, n, n1, n2, k, k1, k2, alpha(ite), th1(ite), th2(ite));

% Test Y
test2a = Yk * [fxy_vec;gxy_vec];
test2b = Sk_fg * x;
test2 = norm(test2a - test2b);
display(test2);

% Calculate the initial residual r = ck - (Ak*x)
res_vec = ck - (Sk_fg*M*xk);

% Get number of entries in the vector f of all perturbations solved by LSE.
nEntries = nNonZeros_fxy...
    + nNonZeros_gxy ...
    + nCoefficients_x ...
    + 3;

% Set the intial value of E to the identity matrix
E = eye(nEntries);

%
% Create the matrix (T+N), initially N is empty so this is the same as T.
Tf = BuildT1_Both_Bivariate(fww, m, n-k, n1-k1, n2-k2);
Tg = BuildT1_Both_Bivariate(gww, n, m-k, m1-k1, m2-k2);
TN = [Tf alpha(ite).*Tg];

%
% Create The matrix (T+N) with respect to alpha
Tfww_wrt_alpha = BuildT1_Both_Bivariate(Partial_fww_wrt_alpha, m, n-k, n1-k1, n2-k2);
Tgww_wrt_alpha = BuildT1_Both_Bivariate(Partial_alpha_gww_wrt_alpha, n, m-k, m1-k1, m2-k2);
TN_wrt_alpha = [Tfww_wrt_alpha Tgww_wrt_alpha];

%
% Create The matrix (T+N) with respect to theta1
Tfww_wrt_th1 = BuildT1_Both_Bivariate(Partial_fww_wrt_th1, m, n-k, n1-k1, n2-k2);
Tgww_wrt_th1 = BuildT1_Both_Bivariate(Partial_gww_wrt_th1, n, m-k, m1-k1, m2-k2);
TN_wrt_theta1 = [Tfww_wrt_th1 alpha(ite) * Tgww_wrt_th1];

%
% Create The matrix (T+N) with respect to theta2
Tfww_wrt_th2 = BuildT1_Both_Bivariate(Partial_fww_wrt_th2, m, n-k, n1-k1, n2-k2);
Tgww_wrt_th2 = BuildT1_Both_Bivariate(Partial_gww_wrt_th2, n, m-k, m1-k1, m2-k2);
TN_wrt_theta2 = [Tfww_wrt_th2 alpha(ite) * Tgww_wrt_th2];

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

yy = start_point;

f = -(yy - start_point);

% Set the termination criterion to a large value. It will be
% over written later.
condition(ite) = norm(res_vec) ./ norm(ck);

while condition(ite) >(SETTINGS.MAX_ERROR_SNTLN) &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN
    %while   ite < max_iterations
    
    
    % Use the QR decomposition to solve the LSE problem
    % min |y-p| subject to Cy=q
    y = LSE(E, f, C, res_vec);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % obtain the small changes
    
    % % 
    % Break down y into its sections
    
    % Get the coefficients corresponding to f and g
    delta_zk = y(1:nNonZeros_fxy + nNonZeros_gxy ,1);
    
    % Remove the zk coefficients from the list of coefficients
    y(1:nNonZeros_fxy + nNonZeros_gxy) = [];
    
    % Get the coefficients corresponding to x
    delta_xk        = y(1:nCoefficients_x,1);
    
    % Remove them from the list of coefficients
    y(1:nCoefficients_x) = [];
    
    % Get the coefficient corresponding to alpha
    delta_alpha     = y(1:1);
    y(1) = [];
    
    % Get the coefficient corresponding to theta1
    delta_theta1    = y(1:1);
    y(1) = [];
    
    % Get the coefficient corresponding to theta2
    delta_theta2    = y(1:1);
    y(1) = [];
    
    % %
    % Update the variables
    
    % Update variables z_{k}, where z_{k} are perturbations in the
    % coefficients of f and g.
    zk = zk + delta_zk;
    
    % Update x_{k}, where x_{k} is the solution vector, containing
    % coefficients u and v.
    xk = xk + delta_xk;
    
    % Update alpha
    alpha(ite) = alpha(ite-1) + delta_alpha;
    
    % Update theta_{1}
    th1(ite) = th1(ite-1) + delta_theta1;
    
    % Update theta_{2}
    th2(ite) = th2(ite-1) + delta_theta2;
    
    % Obtain polynomials in modified bersntein basis a_{i}\theta^{i}
    
    % Obtain new f(w,w) with improved theta1, and theta2
    fww = GetWithThetas(fxy,th1(ite),th2(ite));
    
    % Obtain new g(w,w) with improved theta1 and theta2
    gww = GetWithThetas(gxy,th1(ite),th2(ite));
    
    % Construct the Sylvester subresultant matrix S.
    Tf = BuildT1_Both_Bivariate(fww, m, n-k, n1-k1, n2-k2);
    Tg = BuildT1_Both_Bivariate(gww, n, m-k, m1-k1, m2-k2);
    Sk_fg = [Tf alpha(ite).*Tg];
    
    % Calculate the partial derivatives of fw and gw with respect to alpha
    Partial_fww_wrt_alpha = zeros(m1+1, m2+1);
    Partial_alpha_gww_wrt_alpha  = gww;
    
    % Calculate the partial derivatives of fw and gw with respect to theta1
    % divide the rows by theta1 and multiply by the old power
    
    % Get the partial derivative of f with respect to theta 1
    temp_mat = diag((0:1:m1)./th1(ite));
    Partial_fww_wrt_th1 = temp_mat * fww;
    
    % Get the partial derivative of g with respect to theta1
    temp_mat = diag((0:1:n1)./th1(ite));
    Partial_gww_wrt_th1 = temp_mat * gww;
    
    % Get the partial derivative of f with respect to theta2
    temp_mat = diag((0:1:m2)./th2(ite));
    Partial_fww_wrt_th2 =  fww * temp_mat;
    
    % Get the partial derivative of g with respect too theta2
    temp_mat = diag((0:1:n2)./th2(ite));
    Partial_gww_wrt_th2 =  gww * temp_mat;
    
    
    % Calculate the Partial derivative of T with respect to alpha.
    Tfww_wrt_alpha = BuildT1_Both_Bivariate(Partial_fww_wrt_alpha, m, n-k, n1-k1, n2-k2);
    Tgww_wrt_alpha = BuildT1_Both_Bivariate(Partial_alpha_gww_wrt_alpha, n, m-k, m1-k1, m2-k2);
    Partial_T_wrt_alpha = [Tfww_wrt_alpha Tgww_wrt_alpha];
    
    
    % Calculate the partial derivative of T with respect to theta1
    Tfww_wrt_th1 = BuildT1_Both_Bivariate(Partial_fww_wrt_th1, m, n-k, n1-k1, n2-k2);
    Tgww_wrt_th1 = BuildT1_Both_Bivariate(Partial_gww_wrt_th1, n, m-k, m1-k1, m2-k2);
    Partial_Tfg_wrt_th1 = [Tfww_wrt_th1 alpha(ite)*Tgww_wrt_th1];
    
    % Calculate the partial derivative of T with respect to theta2
    Tfww_wrt_th2 = BuildT1_Both_Bivariate(Partial_fww_wrt_th2, m, n-k, n1-k1, n2-k2);
    Tgww_wrt_th2 = BuildT1_Both_Bivariate(Partial_gww_wrt_th2, n, m-k, m1-k1, m2-k2);
    
    Partial_Tfg_wrt_th2 = [Tfww_wrt_th2 alpha(ite)*Tgww_wrt_th2];
    
    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ck = Sk_fg*e;
    
    % Calculate the derivatives of c_{k} with respect to \alpha and \theta
    Partial_ck_wrt_alpha     = Partial_T_wrt_alpha*e;
    Partial_ck_wrt_th1       = Partial_Tfg_wrt_th1*e;
    Partial_ck_wrt_th2       = Partial_Tfg_wrt_th2*e;
    
    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    z_fxy_vec      = zk(1:nNonZeros_fxy);
    z_gxy_vec      = zk(nNonZeros_fxy + 1 :end);
    
    % Calculate the entries of z_fw and z_gw
    
    z_fxy = GetAsMatrix_Version1([z_fxy_vec ; zeros(nZeros_fxy,1)], m1, m2);
    z_gxy = GetAsMatrix_Version1([z_gxy_vec ; zeros(nZeros_gxy,1)], n1, n2);
    
    % Get z_fw_mat, by multiplying by thetas
    z_fww = GetWithThetas(z_fxy, th1(ite), th2(ite));
    z_gww = GetWithThetas(z_gxy, th1(ite), th2(ite));
    
    % Calculate the derivatives of z_fw and z_gw with repect to alpha.
    Partial_zfww_wrt_alpha = zeros(m1+1,m2+1);
    Partial_zgww_wrt_alpha = z_gww;
    
    % Calculate the derivative of z_fw with respect to theta1.
    temp_mat = diag((0:1:m1)./th1(ite));
    Partial_zfww_wrt_th1 = temp_mat * z_fww;
    
    % Calculate the derivative of z_fw with respect to theta2
    temp_mat = diag((0:1:m2)./th2(ite));
    Partial_zfww_wrt_th2 = z_fww * temp_mat;
    
    % Calculate the derivative of z_gw with respect ot theta1
    temp_mat = diag((0:1:n1)./th2(ite));
    Partial_zgww_wrt_th1 = temp_mat * z_gww;
    
    % Calculate the deriviate of z_gw with respect to theta2
    temp_mat = diag((0:1:n2)./th2(ite));
    Partial_zgww_wrt_th2 = z_gww * temp_mat;
    
    % Build the coefficient Matrix N = [T(z1) T(z2)], of structured perturbations, with
    % same structure as DTQ.
    N1 = BuildT1_Both_Bivariate(z_fww, m, n-k, n1-k1, n2-k2);
    N2 = BuildT1_Both_Bivariate(z_gww, n, m-k, m1-k1, m2-k2);
    N = [N1 alpha(ite).*N2];
    
    % Build the coefficient matrix N with respect to alpha
    T_zf_wrt_alpha = BuildT1_Both_Bivariate(Partial_zfww_wrt_alpha, m, n-k, n1-k1, n2-k2);
    T_zg_wrt_alpha = BuildT1_Both_Bivariate(Partial_zgww_wrt_alpha, n, m-k, m1-k1, m2-k2);
    Partial_N_wrt_alpha = [T_zf_wrt_alpha T_zg_wrt_alpha];
    
    
    % Calculate the derivatives of DNQ with respect to theta
    T_zf_wrt_th1 = BuildT1_Both_Bivariate(Partial_zfww_wrt_th1, m, n-k, n1-k1, n2-k2);
    T_zg_wrt_th1 = BuildT1_Both_Bivariate(Partial_zgww_wrt_th1, n, m-k, m1-k1, m2-k2);
    Partial_N_wrt_th1 = [T_zf_wrt_th1 alpha(ite).*T_zg_wrt_th1];
    
    % Calculate the derivatives of DNQ with respect to theta
    T_zf_wrt_th2 = BuildT1_Both_Bivariate(Partial_zfww_wrt_th2, m, n-k, n1-k1, n2-k2);
    T_zg_wrt_th2 = BuildT1_Both_Bivariate(Partial_zgww_wrt_th2, n, m-k, m1-k1, m2-k2);
    Partial_N_wrt_th2 = [T_zf_wrt_th2 alpha(ite).*T_zg_wrt_th2];
    
    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{k} the column of S_{k} moved to the RHS
    hk = N*e;
    
    % Calculate the derivative of h with respect to alpha
    hk_alpha = Partial_N_wrt_alpha*e;
    
    % Calculate the derivative of h with respect to theta1
    hk_th1 = Partial_N_wrt_th1*e;
    
    % Calculate the derivative of h with respect to theta2
    hk_th2 = Partial_N_wrt_th2*e;
    
    % Build the matrix (T(f,g)+N(z_{f},z_{g}))
    TN_fww = BuildT1_Both_Bivariate(fww + z_fww, m, n-k, n1-k1, n2-k2);
    TN_gww = BuildT1_Both_Bivariate(gww + z_gww, n, m-k, m1-k1, m2-k2);
    TN = [TN_fww alpha(ite).*TN_gww];
    
    % Calculate the paritial derivative of (T+N) with respect to
    % \alpha
    TN1_wrt_alpha = BuildT1_Both_Bivariate(Partial_fww_wrt_alpha + Partial_zfww_wrt_alpha, m, n-k, n1-k1, n2-k2);
    TN2_wrt_alpha = BuildT1_Both_Bivariate(Partial_alpha_gww_wrt_alpha + Partial_zgww_wrt_alpha, n, m-k, m1-k1, m2-k2);
    TN_alpha = [TN1_wrt_alpha TN2_wrt_alpha];
    
    
    % Calculate the paritial derivative of (T+N) with respect to
    % \theta_{1}
    TN1_wrt_th1 = BuildT1_Both_Bivariate(Partial_fww_wrt_th1 + Partial_zfww_wrt_th1, m, n-k, n1-k1, n2-k2);
    TN2_wrt_th1 = BuildT1_Both_Bivariate(Partial_gww_wrt_th1 + Partial_zgww_wrt_th1, n, m-k, m1-k1, m2-k2);
    TN_th1 = [TN1_wrt_th1 alpha(ite).*TN2_wrt_th1];
    
    % Calculate the paritial derivative of (T+N) with respect to
    % \theta_{2}
    TN1_wrt_th2 = BuildT1_Both_Bivariate(Partial_fww_wrt_th2 + Partial_zfww_wrt_th2, m, n-k, n1-k1, n2-k2);
    TN2_wrt_th2 = BuildT1_Both_Bivariate(Partial_gww_wrt_th2 + Partial_zgww_wrt_th2, n, m-k, m1-k1, m2-k2);
    TN_th2 = [TN1_wrt_th2 alpha(ite).*TN2_wrt_th2];
    
    % Update xk
    % Insert a zero into the least squares solution x_ls so that 
    % S_{k,k1,k2} x = c_{k,k1,k2}. Also so that when x is split into two 
    % vectors x1 and x2. S(x1,x2) [f;g] =  ck.
    first_part = xk(1:(idx_col-1));
    second_part = xk(idx_col:end);
    x = [first_part ; 0 ; second_part];

    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    Yk = BuildY_BothDegree_SNTLN(x, m, m1, m2, n, n1, n2, k, k1, k2, alpha(ite), th1(ite), th2(ite));
    
    % Calculate the matrix P where ck = P * [f,g]
    Pk = BuildP_BothDegree_SNTLN(m, m1, m2, n, n1, n2, k, k1, k2, alpha(ite), th1(ite), th2(ite), idx_col);

    % Get residual as a vector
    res_vec = (ck + hk) - TN*M*xk ;
    
    % Create the matrix C. This is made up of five submatrices, HZ, Hx,
    % H_alpha and H_theta1 and H_theta2.
    
    Hz          = Yk - Pk;
    
    Hx          = TN*M;
    
    H_alpha     = TN_alpha*M*xk - (Partial_ck_wrt_alpha + hk_alpha);
    
    H_theta1    = TN_th1*M*xk - (Partial_ck_wrt_th1 + hk_th1);
    
    H_theta2    = TN_th2*M*xk - (Partial_ck_wrt_th2 + hk_th2);
    
    C = [Hz,Hx,H_alpha,H_theta1, H_theta2];  % the matrix C
    
    % Calculate the new right hand vector
    ek = ck + hk;
 
    
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(res_vec) / norm(ek);
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
end

PlotGraphs_SNTLN();

SETTINGS.LOW_RANK_APPROX_REQ_ITE = ite;


%

% Once iterations are complete, assign fx output, gx output, solution X
% output, alpha output and theta output.

% get the vector zk
zPert_f_vec = zk(1:nNonZeros_fxy);
zPert_f_mat = GetAsMatrix_Version1([zPert_f_vec ; zeros(nZeros_fxy,1)],m1,m2);

zPert_g_vec = zk(nNonZeros_fxy+1:end);
zPert_g_mat = GetAsMatrix_Version1([zPert_g_vec ; zeros(nZeros_gxy,1)],n1,n2);

% Set outputs of low rank approximation

fxy_lr = fxy + zPert_f_mat;

gxy_lr = gxy + zPert_g_mat;

first_part = xk(1:(idx_col-1));
second_part = xk(idx_col:end);
x = [first_part ; -1 ; second_part];

vec_vww = x(1:nNonZeros_vxy);
vec_uww = -1.* x(nNonZeros_vxy+1:end);

vww = GetAsMatrix_Version1([vec_vww ; zeros(nZeros_vxy,1)], n1-k1, n2-k2);
uww = GetAsMatrix_Version1([vec_uww ; zeros(nZeros_uxy,1)], m1-k1, m2-k2);

uxy_lr = GetWithoutThetas(uww, th1(ite), th2(ite));
vxy_lr = GetWithoutThetas(vww, th1(ite), th2(ite));

alpha_lr = alpha(ite);

th1_lr = th1(ite);

th2_lr = th2(ite);

    


% Print the number of iterations
LineBreakLarge();
fprintf([mfilename sprintf('SNTLN converged within %i iterations \n', ite)]);
LineBreakLarge();

end


















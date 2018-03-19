function [ fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr,th1_lr,th2_lr] = ...
    SNTLN_Relative( fxy, gxy, i_alpha, i_th1, i_th2, k1, k2, idx_col)
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
% k1 : (Int) Degree of AGCD d(x,y) with respect to x
%
% k2 : (Int) Degree of AGCD d(x,y) with repsect to y
%
% idx_col : (Int) Optimal column for removal from the sylvester matrix, such that col
%           is the column which is most likely a linear combination of the others.
%
%
% Outputs:
%
%
% fxy_lr : (Matrix) Coefficients of f(x,y) on output, in standard bernstein basis,
% including added structured perturbations.
%
% gxy_lr : (Matrix) Coefficients of g(x,y) on output, in standard bernstein basis,
% including added structured perturbations.
%
% uxy_lr : (Matrix) Coefficients of polynomial u(x,y)
%
% vxy_lr : (Matrix) Coefficients of polynomial v(x,y)
%
% alpha_lr : (Float)
%
% th1_lr : (Float) 
%
% th2_lr : (Float) 


% Global Inputs

global SETTINGS

% Set the initial iterations number
ite = 1;

% Set initial values of alpha and theta
th1(ite) = i_th1;
th2(ite) = i_th2;
alpha(ite) = i_alpha;

% Get degree of polynomials f.
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degree of polynomial g.
[n1, n2] = GetDegree_Bivariate(gxy);

% Get the number of coefficients in the polynomial f(x,y)
nCoefficients_fxy = (m1+1) * (m2+1);

% Get the number of coefficients in the polynomial g(x,y)
nCoefficients_gxy = (n1+1) * (n2+1);

% Get the number of coefficients in both f(x,y) and g(x,y)
nCoefficients_fg = nCoefficients_fxy + nCoefficients_gxy;

% Get the number of coefficients in v(x,y)
nCoefficients_vxy = (n1-k1+1) * (n2-k2+1);

% Get the number of coefficients in u(x,y)
nCoefficients_uxy = (m1-k1+1) * (m2-k2+1);

% Get the number of coefficients in the unknown vector x, where A_{t}x =
% c_{t}.
nCoefficients_x = nCoefficients_uxy + nCoefficients_vxy - 1;

% Create the identity matrix I, the matrix M formed from I by removing the
% column equivalent to the optimal column for removal from the Sylvester
% subresultant matrix, so that S_{t}(f,g)*M = A_{t}, where A_{t} is the
% Sylvester subresultant matrix with the column removed.

% Get the number of columns in C_{t}(f), the first partition of the Sylvester
% Matrix S_{t}(f,g)
nColumns_Tf = nCoefficients_vxy;

% Get the number of columns in C_{t}(g), the second partition of the
% Sylvester matrix S_{t}(f,g)
nColumns_Tg = nCoefficients_uxy;

% Get the total number of columns in the Sylvester matrix S_{t}(f,g)
nColumns_Sk = nColumns_Tf + nColumns_Tg;

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

% Preprocessing
% Obtain polynomials in Modified Bernstein Basis, using initial values of
% alpha and theta.

% Multiply the rows of fxy_matrix by theta1, and multiply the cols of
% fxy_matrix by theta2
fww = GetWithThetas(fxy, th1(ite), th2(ite));
gww = GetWithThetas(gxy, th1(ite), th2(ite));

% Form the Coefficient Matrix T = [C(f(w1,w2))|alpha * C(g(w1,w2))] such that T*x = [col]
T1_fww = BuildT1_Relative_Bivariate(fww,n1-k1,n2-k2);
T2_gww = BuildT1_Relative_Bivariate(gww,m1-k1,m2-k2);
T_fg = [T1_fww alpha(ite).*T2_gww];

%
% Calculate the partial derivatives of f(w,w) and g(w,w) with respect to \alpha
Partial_fw_wrt_alpha            = zeros(m1+1,m2+1);
Partial_alpha_gw_wrt_alpha      = gxy;

%
% Calculate the partial derivatives of f(w,w) with respect to \theta_1
theta_mat = diag((0:1:m1) ./ th1(ite));
Partial_fww_wrt_th1    = theta_mat * fww;

%
% Calculate the partial derivative of g(w,w) with respect to theta_1
theta_mat = diag((0:1:n1) ./ th1(ite));
Partial_gw_wrt_th1    = theta_mat * gww;

%
% Calculate the partial derivative of f(w,w) with respect to theta_2
theta_mat = diag((0:1:m2) ./ th2(ite));
Partial_fww_wrt_th2 = fww * theta_mat;

%
% Calculate the partial deriviates of g(w,w) with respect to theta_2
theta_mat = diag((0:1:n2)./ th2(ite));
Partial_gww_wrt_th2 = gww * theta_mat;

%
% Build the derivative of T(f,g) with respect to alpha
T1_fww_wrt_alpha = BuildT1_Relative_Bivariate(Partial_fw_wrt_alpha,n1-k1,n2-k2);
T2_gww_wrt_alpha = BuildT1_Relative_Bivariate(Partial_alpha_gw_wrt_alpha,m1-k1,m2-k2);
T_fg_wrt_alpha = [T1_fww_wrt_alpha T2_gww_wrt_alpha];

%
% Calculate the derivative of T(f,g) with respect to theta_{1}
T1_fww_wrt_th1 = BuildT1_Relative_Bivariate(Partial_fww_wrt_th1,n1-k1,n2-k2);
T2_gww_wrt_th1 = BuildT1_Relative_Bivariate(Partial_gw_wrt_th1,m1-k1,m2-k2);
Partial_T_fg_wrt_th1 = [T1_fww_wrt_th1 alpha(ite)* T2_gww_wrt_th1];

% %
% Calcualte the derivative of T(f,g) with respect to theta_2
T1_fww_wrt_th2 = BuildT1_Relative_Bivariate(Partial_fww_wrt_th2,n1-k1,n2-k2);
T2_gww_wrt_th2 = BuildT1_Relative_Bivariate(Partial_gww_wrt_th2,m1-k1,m2-k2);
Partial_T_fg_wrt_th2 = [T1_fww_wrt_th2 alpha(ite)*T2_gww_wrt_th2];

% %
% Initialise the vector z of structured perturbations
zk = zeros(nCoefficients_fg , 1);

nRows_Sylv_mat = (m1 + n1 - k1 + 1) * (m2 + n2 - k2 + 1);

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
Partial_h_wrt_alpha = Partial_N_wrt_alpha*e;
Partial_h_wrt_th1 = Partial_N_wrt_th1*e;
Partial_h_wrt_th2 = Partial_N_wrt_th2*e;


% Get the matrix A_{k}(f,g), which is the subresultant matrix S(f,g) with
% an opitmal column removed
Ak = T_fg;
ck = T_fg(:,idx_col);
Ak(:,idx_col) = [];

% Build the matrix P
Pk = BuildP_RelativeDegree_SNTLN(m1,m2,n1,n2,k1,k2,alpha(ite),th1(ite),th2(ite),idx_col);
% Test P
% Get the coefficients of f(x,y) in matrix form
fxy_vec = GetAsVector(fxy);
gxy_vec = GetAsVector(gxy);

test1a = Pk*[fxy_vec;gxy_vec];
test1b = ck;
test1 = (test1a - test1b);
display(norm(test1));

%
% Calculate the derivatives of ck wrt alpha and theta.
Partial_ck_wrt_alpha = T_fg_wrt_alpha*e;
Partial_ck_wrt_th1 = Partial_T_fg_wrt_th1*e;
Partial_ck_wrt_th2 = Partial_T_fg_wrt_th2*e;


% %
% Perform QR decomposition of Ak to obtain the solution x
xk = SolveAx_b(Ak,ck);


% % Build Matrix Y, where Y(v,u)*[f;g] = S(f,g)*[u;v]

% Get the vector x such that Sk(f,g)*x = ck
first_part = xk(1:(idx_col-1));
second_part = xk(idx_col:end);
x = [first_part ; 0 ; second_part];

Yk = BuildY_RelativeDegree_SNTLN(x, m1, m2, n1, n2, k1, k2, alpha(ite), th1(ite), th2(ite));

% Test Y
test2a = Yk * [fxy_vec;gxy_vec];
test2b = T_fg * x;
test2 = (test2a - test2b);
display(norm(test2));

% Calculate the initial residual r = ck - (Ak*x)
res_vec = ck - (T_fg*M*xk);

% Get the matrix p, which will store all the perturbations returned from LSE file
nEntries = nCoefficients_fxy...
    + nCoefficients_gxy ...
    + nCoefficients_x ...
    + 3;

%
% Set the intial value of E to the identity matrix
E = eye(nEntries);

%
% Create the matrix (T+N), initially N is empty so this is the same as T.
T1_fww = BuildT1_Relative_Bivariate(fww,n1-k1,n2-k2);
T2_gww = BuildT1_Relative_Bivariate(gww,m1-k1,m2-k2);
TN = [T1_fww alpha(ite).*T2_gww];

%
% Create The matrix (T+N) with respect to \alpha
T1_fww_wrt_alpha = BuildT1_Relative_Bivariate(Partial_fw_wrt_alpha,n1-k1,n2-k2);
T2_gww_wrt_alpha = BuildT1_Relative_Bivariate(Partial_alpha_gw_wrt_alpha,m1-k1,m2-k2);
TN_wrt_alpha = [T1_fww_wrt_alpha T2_gww_wrt_alpha];

%
% Create The matrix (T+N) with respect to \theta_{1}
T1_fww_wrt_th1 = BuildT1_Relative_Bivariate(Partial_fww_wrt_th1,n1-k1,n2-k2);
T2_gww_wrt_th1 = BuildT1_Relative_Bivariate(Partial_gw_wrt_th1,m1-k1,m2-k2);
TN_wrt_theta1 = [T1_fww_wrt_th1 alpha(ite) * T2_gww_wrt_th1];

%
% Create The matrix (T+N) with respect to \theta_{2}
T1_fww_wrt_th2 = BuildT1_Relative_Bivariate(Partial_fww_wrt_th2,n1-k1,n2-k2);
T2_gww_wrt_th2 = BuildT1_Relative_Bivariate(Partial_gww_wrt_th2,m1-k1,m2-k2);
TN_wrt_theta2 = [T1_fww_wrt_th2 alpha(ite) * T2_gww_wrt_th2];

%
% Create the matrix C for input into iteration

H_z     = Yk-Pk;

H_x     = TN*M;

H_alpha  = TN_wrt_alpha*M*xk - ...
    (Partial_ck_wrt_alpha + Partial_h_wrt_alpha);

H_theta1 = TN_wrt_theta1*M*xk - ...
    (Partial_ck_wrt_th1 + Partial_h_wrt_th1);

H_theta2 = TN_wrt_theta2*M*xk - ...
    (Partial_ck_wrt_th2 + Partial_h_wrt_th2);

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

f = -(yy-start_point);

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
    
    % Break down y into its sections
    
    % Get the coefficients corresponding to f and g
    delta_zk        = y(1:nCoefficients_fxy + nCoefficients_gxy ,1);
    
    % Remove the zk coefficients from the list of coefficients
    y(1:nCoefficients_fxy + nCoefficients_gxy) = [];
    
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
    
    % Obtain new f(w,w) and g(w,w)
    fww = GetWithThetas(fxy,th1(ite),th2(ite));
    gww = GetWithThetas(gxy,th1(ite),th2(ite));
    
    % Construct the Sylvester subresultant matrix S.
    T1_fww = BuildT1_Relative_Bivariate(fww,n1-k1,n2-k2);
    T2_gww = BuildT1_Relative_Bivariate(gww,m1-k1,m2-k2);
    T_fg = [T1_fww alpha(ite).*T2_gww];
    
    % Calculate the partial derivatives of fw and gw with respect to alpha
    Partial_fw_wrt_alpha            = zeros(m1+1,m2+1);
    Partial_alpha_gw_wrt_alpha      = gww;
    
    %
    % Calculate the partial derivatives of fw and gw with respect to theta1
    % divide the rows by theta1 and multiply by the old power
    
    % Get the partial derivative of f with respect to \theta_{1}
    temp_mat = diag((0:1:m1)./th1(ite));
    Partial_fww_wrt_th1 = temp_mat * fww;
    
    % Get the partial derivative of g with respect to \theta_{1}
    temp_mat = diag((0:1:n1)./th1(ite));
    Partial_gww_wrt_th1 = temp_mat * gww;
    
    % Get the partial derivative of f with respect to \theta_{2}
    temp_mat = diag((0:1:m2)./th2(ite));
    Partial_fww_wrt_th2 =  fww * temp_mat;
    
    % Get the partial derivative of g with respect to \theta_{2}
    temp_mat = diag((0:1:n2)./th2(ite));
    Partial_gww_wrt_th2 =  gww * temp_mat;
    
    
    % Calculate the Partial derivative of T with respect to \alpha.
    T1_fww_wrt_alpha = BuildT1_Relative_Bivariate(Partial_fw_wrt_alpha, n1-k1, n2-k2);
    T2_gww_wrt_alpha = BuildT1_Relative_Bivariate(Partial_alpha_gw_wrt_alpha, m1-k1, m2-k2);
    Partial_T_fg_wrt_alpha = [T1_fww_wrt_alpha T2_gww_wrt_alpha];
    
    % Calculate the partial derivative of T with respect to \theta_{1}
    T1_fww_wrt_th1 = BuildT1_Relative_Bivariate(Partial_fww_wrt_th1, n1-k1, n2-k2);
    T2_gww_wrt_th1 = BuildT1_Relative_Bivariate(Partial_gww_wrt_th1, m1-k1, m2-k2);
    Partial_T_fg_wrt_th1 = [T1_fww_wrt_th1 alpha(ite)*T2_gww_wrt_th1];
    
    % Calculate the partial derivative of T with respect to \theta_{2}
    T1_fww_wrt_th2 = BuildT1_Relative_Bivariate(Partial_fww_wrt_th2, n1-k1, n2-k2);
    T2_gww_wrt_th2 = BuildT1_Relative_Bivariate(Partial_gww_wrt_th2, m1-k1, m2-k2);
    
    Partial_T_fg_wrt_th2 = [T1_fww_wrt_th2 alpha(ite)*T2_gww_wrt_th2];
    
    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ck = T_fg*e;
    
    % Calculate the derivatives of c_{k} with respect to \alpha and \theta
    Partial_ck_wrt_alpha = Partial_T_fg_wrt_alpha*e;
    Partial_ck_wrt_th1 = Partial_T_fg_wrt_th1*e;
    Partial_ck_wrt_th2 = Partial_T_fg_wrt_th2*e;
    
    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    z_fxy      = zk(1:nCoefficients_fxy);
    z_gxy      = zk(nCoefficients_fxy + 1 :end);
    
    % Calculate the entries of z_fw and z_gw
    z_fxy = GetAsMatrix(z_fxy,m1,m2);
    z_gxy = GetAsMatrix(z_gxy,n1,n2);
    
    % Get z_fw_mat, by multiplying by thetas
    z_fww = GetWithThetas(z_fxy,th1(ite),th2(ite));
    z_gww = GetWithThetas(z_gxy,th1(ite),th2(ite));
    
    % Calculate the derivatives of z_fw and z_gw with repect to alpha.
    Partial_zfww_wrt_alpha    = zeros(m1+1,m2+1);
    Partial_zgww_wrt_alpha    = z_gww;
    
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
    N1 = BuildT1_Relative_Bivariate(z_fww,n1-k1,n2-k2);
    N2 = BuildT1_Relative_Bivariate(z_gww,m1-k1,m2-k2);
    N = [N1 alpha(ite).*N2];
    
    % Build the coefficient matrix N with respect to alpha
    T1_zf_wrt_alpha = BuildT1_Relative_Bivariate(Partial_zfww_wrt_alpha, n1-k1,n2-k2);
    T2_zg_wrt_alpha = BuildT1_Relative_Bivariate(Partial_zgww_wrt_alpha, m1-k1,m2-k2);
    Partial_N_wrt_alpha = [T1_zf_wrt_alpha T2_zg_wrt_alpha];
    
    % Calculate the derivatives of DNQ with respect to theta
    T1_zf_wrt_th1 = BuildT1_Relative_Bivariate(Partial_zfww_wrt_th1,n1-k1,n2-k2);
    T2_zg_wrt_th1 = BuildT1_Relative_Bivariate(Partial_zgww_wrt_th1,m1-k1,m2-k2);
    Partial_N_wrt_th1 = [T1_zf_wrt_th1 alpha(ite).*T2_zg_wrt_th1];
    
    % Calculate the derivatives of DNQ with respect to theta
    T1_zf_wrt_th2 = BuildT1_Relative_Bivariate(Partial_zfww_wrt_th2,n1-k1,n2-k2);
    T2_zg_wrt_th2 = BuildT1_Relative_Bivariate(Partial_zgww_wrt_th2,m1-k1,m2-k2);
    Partial_N_wrt_th2 = [T1_zf_wrt_th2 alpha(ite).*T2_zg_wrt_th2];
    
    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{k} the column of S_{k} moved to the RHS
    hk = N*e;
    
    % Calculate the derivative of h_{k} with respect to \alpha
    Partial_hk_alpha = Partial_N_wrt_alpha*e;
    
    % Calculate the derivative of h with respect to \theta_{1}
    Partial_hk_th1 = Partial_N_wrt_th1*e;
    
    % Calculate the derivative of h with respect to \theta_{2}
    Partial_hk_th2 = Partial_N_wrt_th2*e;
    
    % Build the matrix (T+N)
    
    TN1_fww = BuildT1_Relative_Bivariate(fww + z_fww,n1-k1,n2-k2);
    TN2_gww = BuildT1_Relative_Bivariate(gww + z_gww,m1-k1,m2-k2);
    TN = [TN1_fww alpha(ite).*TN2_gww];
    
    % Calculate the paritial derivative of (T+N) with respect to
    % alpha
    TN1_wrt_alpha = BuildT1_Relative_Bivariate(Partial_fw_wrt_alpha + Partial_zfww_wrt_alpha, n1-k1,n2-k2);
    TN2_wrt_alpha = BuildT1_Relative_Bivariate(Partial_alpha_gw_wrt_alpha + Partial_zgww_wrt_alpha, m1-k1,m2-k2);
    TN_alpha = [TN1_wrt_alpha TN2_wrt_alpha];
    
    
    % Calculate the paritial derivative of (T+N) with respect to
    % theta1
    TN1_wrt_th1 = BuildT1_Relative_Bivariate(Partial_fww_wrt_th1 + Partial_zfww_wrt_th1,n1-k1,n2-k2);
    TN2_wrt_th1 = BuildT1_Relative_Bivariate(Partial_gww_wrt_th1 + Partial_zgww_wrt_th1,m1-k1,m2-k2);
    TN_th1 = [TN1_wrt_th1 alpha(ite).*TN2_wrt_th1];
    
    % Calculate the paritial derivative of (T+N) with respect to
    % theta2
    TN1_wrt_th2 = BuildT1_Relative_Bivariate(Partial_fww_wrt_th2 + Partial_zfww_wrt_th2,n1-k1,n2-k2);
    TN2_wrt_th2 = BuildT1_Relative_Bivariate(Partial_gww_wrt_th2 + Partial_zgww_wrt_th2,m1-k1,m2-k2);
    TN_th2 = [TN1_wrt_th2 alpha(ite).*TN2_wrt_th2];
    
   
    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.    
    first_part = xk(1:(idx_col-1));
    second_part = xk(idx_col:end);
    x = [first_part ; 0 ; second_part];

    Yk = BuildY_RelativeDegree_SNTLN(x, m1, m2, n1, n2, k1, k2, alpha(ite), th1(ite), th2(ite));
    
    % Test Y
    
    % Calculate the matrix P where ck = P * [f,g]
    Pk = BuildP_RelativeDegree_SNTLN(m1, m2, n1, n2, k1, k2, alpha(ite), th1(ite), th2(ite), idx_col);
    
    % Get residual as a vector
    res_vec = (ck + hk) - (TN*M*xk) ;
    
    % Create the matrix C. This is made up of five submatrices, HZ, Hx,
    % H_alpha and H_theta1 and H_theta2.
    
    Hz          = Yk - Pk;
    
    Hx          = TN*M;
    
    H_alpha     = TN_alpha*M*xk - (Partial_ck_wrt_alpha + Partial_hk_alpha);
    
    H_theta1    = TN_th1*M*xk - (Partial_ck_wrt_th1 + Partial_hk_th1);
    
    H_theta2    = TN_th2*M*xk - (Partial_ck_wrt_th2 + Partial_hk_th2);
    
    C = [Hz, Hx, H_alpha, H_theta1, H_theta2];  % the matrix C
    
    % Calculate the new right hand vector
    ek = ck + hk;
    
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(res_vec) / norm(ek);
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
end

PlotGraphs_SNTLN()
SETTINGS.LOW_RANK_APPROX_REQ_ITE = ite;

% Once iterations are complete, assign fx output, gx output, solution X
% output, alpha output and theta output.

% Get the vector zk
zPert_f_vec = zk(1:nCoefficients_fxy);
zPert_f_mat = GetAsMatrix(zPert_f_vec,m1,m2);

zPert_g_vec = zk(nCoefficients_fxy+1:end);
zPert_g_mat = GetAsMatrix(zPert_g_vec,n1,n2);

% Set outputs of low rank approximation

fxy_lr = fxy + zPert_f_mat;

gxy_lr = gxy + zPert_g_mat;

% %
% Get coefficients of u(x,y) and v(x,y)
first_part = xk(1:(idx_col-1));
second_part = xk(idx_col:end);
x = [first_part ; -1 ; second_part];
vec_vww = x(1:nCoefficients_vxy);
vec_uww = -1.* x(nCoefficients_vxy+1:end);

% Get matrix of coefficients of v(\omega_{1},\omega_{2})
vww = GetAsMatrix(vec_vww,n1-k1,n2-k2);

% Get matrix of coefficients of u(\omega_{1},\omega_{2}0
uww = GetAsMatrix(vec_uww,m1-k1,m2-k2);

% Get u(x,y) and v(x,y)
vxy_lr = GetWithoutThetas(vww,th1(ite),th2(ite));
uxy_lr = GetWithoutThetas(uww,th1(ite),th2(ite));

% %
% Get \alpha, \theta_{1} and \theta_{2}
alpha_lr = alpha(ite);

% Get \theta_{1}
th1_lr = th1(ite);

% Get \theta_{2}
th2_lr = th2(ite);


% Print the number of iterations
LineBreakLarge();
fprintf([mfilename ' : ' sprintf('SNTLN converged within %i iterations \n', ite)]);
LineBreakLarge()

end


















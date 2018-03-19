function [fx_lr, gx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = ...
    SNTLN(fx, gx, alpha, theta, k)
% Given two input polynomials and the degree of their GCD, Obtain the Low
% Rank Approximation Sylvester Matrix

% Globals
global SETTINGS


% Get degree of polynomial f(x)
m = GetDegree(fx);

% Get degree of polynomial g(x)
n = GetDegree(gx);

% Set the initial iteration number
ite = 1;

% Create the identity matrix I, such that S*I = S
I = eye(m+n-2*k+2);

% Get f(w) from f(x)
fw = GetWithThetas(fx,theta);

% Get partial f(w) wrt \theta
fw_wrt_theta = Differentiate_wrt_theta(fw,theta);

% Get partial f(w) wrt \alpha
fw_wrt_alpha = zeros(m+1,1);

% Get g(w)
gw = GetWithThetas(gx,theta);

% Get partial g(w) wrt \theta
gw_wrt_theta = Differentiate_wrt_theta(gw,theta);

% Get partial g(w) wrt \alpha
gw_wrt_alpha = gw;

% Get the matrix T, used to construct S_{t}, where S_{t} = DTQ.
T = BuildT(fw,alpha.*gw,k);

% Get the index of the optimal colummn for removal
[~,idx_col] = GetMinDistance(T);

% Create the matrix M, such that S_{t}*M = A_{t}
M = I;
M(:,idx_col) = [];

% Create the matrix e, such that S_{t}*e = c_{t}
e = I(:,idx_col);

% Get the partial T wrt alpha
T_wrt_alpha = BuildT(fw_wrt_alpha,gw_wrt_alpha,k);

% Get partial T wrt theta
T_wrt_theta = BuildT(fw_wrt_theta, alpha.*gw_wrt_theta,k);

% Get the Matrix A_{t} = Sk with column c_{t} removed
At = T*M;

% Get the column c_{t}, the removed column from S_{t}
ck = T*e;

% Get the column partial c_{t} wrt \alpha
ck_wrt_alpha = T_wrt_alpha(:,idx_col);

% Get the column partial c_{t} wrt \theta
ck_wrt_theta = T_wrt_theta(:,idx_col);

% Initialise the vector of perturbations z, corresponding to polynomials
% f and g.
z = zeros(m+n+2,1);
z_fx = zeros(m+1,1);
z_gx = zeros(n+1,1);

% Build the matrix N where N has the same structure as T
N = BuildT(z_fx,alpha.*z_gx,k);

% Get the partial derivative of N_{t}
N_wrt_theta = zeros(m+n-k+1, m+n-(2*k)+2);
N_wrt_alpha = zeros(m+n-k+1, m+n-(2*k)+2);

% Get the column h_{t} - s.t h_{t} has the same structure as c_{t}
hk = N*e;
h_wrt_alpha     = N_wrt_alpha*e;
h_wrt_theta     = N_wrt_theta*e;

% Build the matrix T + N
TN = T + N;

TN_wrt_alpha = T_wrt_alpha + N_wrt_alpha;
TN_wrt_theta = T_wrt_theta + N_wrt_theta;

% Calculate the initial estimate of x - the vector whcih contains the
% coefficients of the quotient polynomials u and v.
xk = SolveAx_b(At,ck);

x_a = xk(1:idx_col-1);
x_b = xk(idx_col:end);
% Get the vector x(w), where x includes thetas
x = [x_a ;0 ;x_b];

% Build the matrix P_{t} where P_{t}*[f;g] = c_{t}
Pk = BuildP_SNTLN(m,n,k,alpha(ite),theta(ite),idx_col);

% Build the matrix Y such that Y*[f;g] = A*x
Yk = BuildY_SNTLN(m,n,k,x,alpha(ite),theta(ite));

% Create the matrix C for input into iteration

H_z = Yk - Pk;

H_x = TN*M;

H_alpha = TN_wrt_alpha*M*xk - ...
    (ck_wrt_alpha + h_wrt_alpha);

H_theta = TN_wrt_theta*M*xk - ...
    (ck_wrt_theta + h_wrt_theta);

C = [H_z H_x H_alpha H_theta];

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    z;...
    xk;...
    alpha(ite);...
    theta(ite);
    ];

yy  = start_point;

%%
% Create matrix E.
E = eye(2*m+2*n-2*k+5);

% Set the initial value of vector p to be zero
% f = -(yy - start_point);
f = zeros(2*m+2*n-2*k+5,1);

% Get the initial residual vector
res_vec = (ck+hk) - (At*xk);

% Get residual
condition = norm(res_vec) ./ norm(ck);


while condition(ite) > SETTINGS.MAX_ERROR_SNTLN && ite < SETTINGS.MAX_ITE_SNTLN
    
    % Perfome LSE Calculation min|Ey-f| subject to Cy=g
    y_lse = LSE(E,f,C,res_vec);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y_lse;
    
    % obtain the small changes
    delta_zk        = y_lse(1:m+n+2,1);
    delta_xk        = y_lse((m+n+3):(2*m+2*n-2*k+3),1);
    delta_alpha     = y_lse(2*m+2*n-2*k+4);
    delta_theta     = y_lse(2*m+2*n-2*k+5);
    
    % Update variables z_{k}, x_{k}, where z_{k} are perturbations in the
    % coefficients of f and g. x_{k} is the solution vector, containing
    % coefficients u and v.
    z = z + delta_zk;
    
    xk = xk + delta_xk;

    % Update \alpha and \theta
    alpha(ite) = alpha(ite-1) + delta_alpha;
    theta(ite) = theta(ite-1) + delta_theta;
    
    % Obtain polynomials in modified basis f(w,\theta)
    fw = GetWithThetas(fx,theta(ite));
    gw = GetWithThetas(gx,theta(ite));
    
    % Construct the subresultant matrix of T.
    T = BuildT(fw,alpha(ite).*gw,k);
    
    % Get new c_{t}
    ck = T*e;
    
    % Calculate the partial derivatives of fw and gw with respect to alpha
    fw_wrt_alpha    = zeros(m+1,1);
    gw_wrt_alpha    = gw;
    
    % Calculate the partial derivatives of fw and gw with respect to theta
    fw_wrt_theta = Differentiate_wrt_theta(fw,theta(ite));
    gw_wrt_theta = Differentiate_wrt_theta(gw,theta(ite));
    
    % Calculate the partial derivative of T wrt alpha
    T_wrt_alpha = BuildT(fw_wrt_alpha,gw_wrt_alpha,k);
    
    % Calculate the partial derivative of T wrt theta
    T_wrt_theta = BuildT(fw_wrt_theta, alpha(ite).*gw_wrt_theta,k);
    
    
    % Calculate the derivatives of c_{k} with respect to \alpha and \theta
    ck_wrt_alpha     = T_wrt_alpha*e;
    ck_wrt_theta     = T_wrt_theta*e;
    
    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    z_fx      = z(1:m+1);
    z_gx      = z(m+2:end);
    
    % Calculate the subresultant matrix of the structured perturbations.
    z_fw = GetWithThetas(z_fx,theta(ite));
    z_gw = GetWithThetas(z_gx,theta(ite));
    
    % Calculate the derivatives of z_fw and z_gw with repect to alpha.
    partial_zfw_wrt_alpha    = zeros(m+1,1);
    partial_zgw_wrt_alpha    = z_gw;
    
    % Calculate the derivatives of z_fw and z_gw with respect to theta.
    partial_zfw_wrt_theta = Differentiate_wrt_theta(z_fw,theta(ite));
    partial_zgw_wrt_theta = Differentiate_wrt_theta(z_gw,theta(ite));
    
    % Build the Coefficient Matrix N, of structured perturbations, with
    % same structure as T.
    N = BuildT(z_fw,alpha(ite).*z_gw,k);
    
    % Build the Sylvester matrix
    N_wrt_alpha = BuildT(partial_zfw_wrt_alpha, partial_zgw_wrt_alpha,k);
    
    % Build the Syvlester matrix
    N_wrt_theta = BuildT(partial_zfw_wrt_theta, alpha(ite).*partial_zgw_wrt_theta,k);
    
    % update h - the vector of structured perturbations equivalent to ck
    hk = N*e;
    
    % Calculate the derivative of h with respect to alpha
    hk_alpha = N_wrt_alpha*e;
    
    % Calculate the derivative of h with respect to theta
    hk_theta = N_wrt_theta*e;
    
    % Update T+N
    TN = T + N;
    
    % Update parital derivative of T+N wrt alpha
    TN_alpha = T_wrt_alpha + N_wrt_alpha;
    
    % Update parital derivative of T+N wrt theta
    TN_theta = T_wrt_theta + N_wrt_theta;
    
    
    x_a = xk(1:idx_col-1);
    x_b = xk(idx_col:end);
    
    % Get the vector x(w), where x includes thetas
    x = [x_a ;0 ;x_b];
    
    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    Yk = BuildY_SNTLN(m,n,k,x,alpha(ite),theta(ite));
    
    % Calculate the matrix P where P is the matrix such that c = P[f;g]
    Pk = BuildP_SNTLN(m,n,k,alpha(ite),theta(ite),idx_col);
    
    % Calculate the residual g
    res_vec = (ck+hk) - ((TN * M) * xk);
    
    % Create the matrix C. This is made up of four submatrices, HZ, Hx,
    % H_alpha and H_theta.
    
    Hz = Yk - Pk;
    
    Hx = TN*M;
    
    H_alpha = TN_alpha*M*xk - (ck_wrt_alpha + hk_alpha);
    
    H_theta = TN_theta*M*xk - (ck_wrt_theta + hk_theta);
    
    C = [Hz,Hx,H_alpha,H_theta];  % the matrix C
    
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(res_vec) ./ norm(ck+hk) ;
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    
end

LineBreakLarge()
fprintf('Required Number of iterations : %i \n',ite)
LineBreakLarge()
SETTINGS.LOW_RANK_APPROX_REQ_ITE = ite;

fx_lr = fx + z_fx;
gx_lr = gx + z_gx;

x_a = xk(1:idx_col-1);
x_b = xk(idx_col:end);

% Get the vector x(w), where x includes thetas
x = [x_a ; -1 ;x_b];

% Get number of coefficients in the polynomial v(x,y)
nCoeffs_vx = n-k+1;

% Get coefficients of v(\omega)
vw = x(1:nCoeffs_vx);

% Get coefficients of u(\omega)
uw = -1.*x(nCoeffs_vx + 1 : end);

% Get coefficients of u(x) and v(x)
ux_lr = GetWithoutThetas(uw,theta(ite));
vx_lr = GetWithoutThetas(vw,theta(ite));

% Get \theta and \alpha
theta_lr = theta(ite);
alpha_lr = alpha(ite);


end






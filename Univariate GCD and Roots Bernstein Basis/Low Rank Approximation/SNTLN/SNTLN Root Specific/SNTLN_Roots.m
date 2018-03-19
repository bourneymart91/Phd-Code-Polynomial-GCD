function [ fx_output,gx_output,alpha_output,theta_output,X_output] = ...
    SNTLN_Roots(fx_n, gx_n, i_alpha, theta, k, idx_col, gm_fx, gm_gx)
%                             Inputs.
%
% fx_n : Coefficients of polynomial f, in standard bernstein basis. where
%        fx is noisy,  and has been normalised by geometric mean
%
% gx_n : Coefficients of polynomial g, in standard bernstein basis. where
%        gx is noisy, and has been normalised by geometric mean
%
%
% opt_col : Optimal column for removal from the sylvester matrix, such that col
%           is the column which is most likely a linear combination of the others.
%
% i_alpha : Initial value of alpha
%
% i_th : Initial value of theta
%
% t  : Degree of AGCD.
%
% gm_f : Geometric mean of entries of f in the Sylvester matrix.
%
% gm_g : Geometric mean of entries of g in the Syvlester matrix.
%
% % Outputs.
%
%
% fx_output :- Coefficients of fx on output, in standard bernstein basis,
% including added structured perturbations.
%
% gx_output :- Coefficients of fx on output, in standard bernstein basis,
% including added structured perturbations.
%
% alpha_output :-
%
% theta_output :-
%
% X_output :-
%
%

% % Global Inputs
% These values should only be reset in the head function o_roots
global SETTINGS

% Get ratio of geometric means.
ratio = gm_fx./gm_gx;

% Set the initial iterations number
ite = 1;

% Set initial values of alpha
alpha(ite) = i_alpha;

% Set initial values of theta
th(ite) = theta;

% Get degree of polynomials f.
m = GetDegree(fx_n);

% Get degree of polynomial g.
n = GetDegree(gx_n);

% Create the identity matrix I
I = eye(m+n-2*k+2,m+n-2*k+2);

% Form matrix M from the identity matrix such that S(f,g)*M gives A_{t},
% the sylvester matrix with the optimal column removed.
M = I;
M(:,idx_col) = [];

% Form vector e such that S(f,g)*e gives the column c_{t}
e = I(:,idx_col);

% Replace the gx which has been input, obtain g(x) as derivative of
% f(x)
gx_new = zeros(m,1);
for i = 1:1:length(fx_n)-1
    gx_new(i) = m.*ratio.*(fx_n(i+1)-fx_n(i));
end

% [gx_new gx_n gx_new-gx_n]

% Obtain polynomials in Modified Bernstein Basis, using initial values of
% alpha and theta.
fw = GetWithThetas(fx_n,th(ite));
gw = GetWithThetas(gx_new,th(ite));

% Form the Coefficient Matrix DTQ such that DTQ * x = [col]
DTQ = BuildDTQ(fw,alpha(ite).*gw,k);

% Calculate the partial derivatives of fw and gw with respect to alpha
fw_wrt_alpha            = zeros(m+1,1);
alpha_gw_wrt_alpha      = gw;

% Calculate the partial derivatives of fw and gw with respect to theta
fw_wrt_theta    = Differentiate_wrt_theta(fw,th(ite));
gw_wrt_theta    = Differentiate_wrt_theta(gw,th(ite));

% Calculate derivative of D_{k}T(f,g)Q_{k} with respect to alpha
DTQ_alpha = BuildDTQ(fw_wrt_alpha,alpha_gw_wrt_alpha,k);

% Calculate the derivative of D_{k}TQ_{k} with respect to theta
DTQ_theta = BuildDTQ(fw_wrt_theta,alpha(ite).*gw_wrt_theta,k);

% Initialise the vector z of structured perturbations
% When working with the roots problem, the number of entries in z is equal
% The number of coefficients of f(x), The perturbations of g(x) are the
% derivative of those from f(x).
zk   = zeros(m+1,1);
z_fx = zeros(m+1,1);
z_gx = zeros(m,1);

% Initilaise the derivative of D_{k}NQ_{k} wrt alpha.
DNQ_wrt_alpha   = zeros(m+n-k+1,m+n-(2*k)+2);

% Initilaise the derivative of D_{k}NQ_{k} wrt theta.
DNQ_wrt_theta   = zeros(m+n-k+1,m+n-(2*k)+2);

% Calculate the derivatives wrt alpha and theta of the column of DNQ
% that is moved to the right hand side.
Partial_h_wrt_alpha     = DNQ_wrt_alpha*e;

Partial_h_wrt_theta     = DNQ_wrt_theta*e;

%Calculate the matrix P. such that Pf = ck
DP  = BuildDP_Roots(m,n,alpha(ite),th(ite),idx_col,k,ratio);

% Get the column c_{t}
ct = (DTQ)*e;

% Calculate the remaining columns of the matrix.
At = (DTQ)*M;

% Calculate the derivatives wrt alpha and theta of the removed column.
ct_wrt_alpha     = DTQ_alpha*e;
ct_wrt_theta     = DTQ_theta*e;

% Calculate the initial estimate of x - the vector whcih contains the
% coefficients of the quotient polynomials u and v.
x_ls = SolveAx_b(At,ct);


% DY0 =   BuildDY(m,n,t,min_col,x_ls,alpha(ite),theta(ite));
DY  =   BuildDY_Roots(m,n,k,idx_col,x_ls,alpha(ite),th(ite),ratio);

% Calculate the initial residual r = ck - (Ak*x)
res_vec = ct - (DTQ*M*x_ls);

% Set the initial value of vector p to be zero

p = zeros(2*m+n-2*k+4,1);

% Set the intial value of E to the identity matrix
E = eye(2*m+2*n-2*k+5);

% Create the matrix D(T+N)Q
DTNQ = BuildDTQ(fw,alpha(ite).*gw,k)

% Create The matrix D(T+N)Q with respect to alpha
DTNQ_wrt_alpha = BuildDTQ(fw_wrt_alpha, alpha_gw_wrt_alpha,k);

% Create The matrix D(T+N)Q with respect to theta
DTNQ_wrt_theta = BuildDTQ(fw_wrt_theta, alpha(ite).*gw_wrt_theta,k);


% Create the matrix C for input into iteration

H_z     = GetH_z(idx_col,n,k,DY,DP);

H_x     = DTNQ*M;

H_beta  = DTNQ_wrt_alpha*M*x_ls - ...
    (ct_wrt_alpha + Partial_h_wrt_alpha);

H_omega = DTNQ_wrt_theta*M*x_ls - ...
    (ct_wrt_theta + Partial_h_wrt_theta);

C       = [H_z H_x H_beta H_omega];

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    zk;...
    x_ls;...
    alpha(ite);...
    th(ite)
    ];

yy      =   start_point;

% Set the termination criterion to a large value. It will be
% over written later.
condition(ite) = norm(res_vec)/norm(ct);

xk = x_ls;

while condition(ite) >(SETTINGS.MAX_ERROR_SNTLN) &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN
    
    % Use the QR decomposition to solve the LSE problem
    % min |y-p| subject to Cy=q
    y = LSE(E,p,C,res_vec);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % obtain the small chages
    delta_zk        = y(1:m+1,1);
    delta_xk        = y((m+2):(2*m+n-2*k+2),1);
    delta_alpha     = y(2*m+n-2*k+3);
    delta_theta     = y(2*m+n-2*k+4);
    
    % Update variables zk, xk
    zk = zk + delta_zk;
    xk = xk + delta_xk;
    
    % Update alpha and theta
    alpha(ite)  = alpha(ite-1) + delta_alpha;
    th(ite)     = th(ite-1) + delta_theta;
    
    % Obtain polynomials in modified bersntein basis a_{i}\theta^{i}
    fw = GetWithThetas(fx_n,th(ite));
    gw = GetWithThetas(gx_new,th(ite));
    
    % Construct the subresultant matrix of DTQ.
    DTQ = BuildDTQ(fw,alpha(ite).*gw,k);
    
    % Calculate the partial derivatives of fw and gw with respect to alpha
    fw_wrt_alpha    = zeros(m+1,1);
    Partial_gw_wrt_alpha    = gw;
    
    % Calculate the partial derivatives of fw and gw with respect to theta
    fw_wrt_theta    = Differentiate_wrt_theta(fw,th(ite));
    gw_wrt_theta    = Differentiate_wrt_theta(gw,th(ite));
    
    % Calculate the Partial derivative of DTQ with respect to alpha.
    Partial_DTQ_wrt_alpha = BuildDTQ(fw_wrt_alpha, Partial_gw_wrt_alpha,k);
    
    % Calculate the partial derivative of DTQ with respect to theta.
    Partial_DTQ_wrt_theta = BuildDTQ(fw_wrt_theta, alpha(ite).*gw_wrt_theta,k);
    
    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ct = DTQ*e;
    
    % Calculate the derivatives of c_{k} with respect to \alpha and \theta
    ct_wrt_alpha     = Partial_DTQ_wrt_alpha*e;
    ct_wrt_theta     = Partial_DTQ_wrt_theta*e;
    
    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    z_fx   = zk(1:m+1);
    
    % EDIT 15/05/2015
    % Get the small changes to gx from the small changes in fx
    z_gx = zeros(m,1);
    for i = 1:1:length(z_fx)-1
        z_gx(i) = m.*ratio.*(z_fx(i+1) - z_fx(i))  ;
    end
    
    % Calculate the subresultant matrix of the structured perturbations.
    z_fw = GetWithThetas(z_fx,th(ite));
    z_gw = GetWithThetas(z_gx,th(ite));
    
    % Calculate the derivatives of z_fw and z_gw with repect to alpha.
    Partial_zfw_wrt_alpha    = zeros(m+1,1);
    Partial_zgw_wrt_alpha    = z_gw;
    
    % Calculate the derivatives of z_fw and z_gw with respect to theta.
    Partial_zfw_wrt_theta    = Differentiate_wrt_theta(z_fw,th(ite));
    Partial_zgw_wrt_theta    = Differentiate_wrt_theta(z_gw,th(ite));
    
    % Build the Coefficient Matrix DNQ, of structured perturbations, with
    % same structure as DTQ.
    DNQ = BuildDTQ(z_fw,alpha(ite).*z_gw,k);
    
    % Calculate the derivatives of DNQ with respect to alpha
    DNQ_wrt_alpha = BuildDTQ(Partial_zfw_wrt_alpha, Partial_zgw_wrt_alpha,k);
    
    % Calculate the derivatives of DNQ with respect to theta
    DNQ_wrt_theta = BuildDTQ(Partial_zfw_wrt_theta, alpha(ite).*Partial_zgw_wrt_theta,k);
    
    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{k} the column of S_{k} moved to the RHS
    h = DNQ*e;
    
    % Calculate the derivative of h with respect to alpha
    h_alpha = DNQ_wrt_alpha*e;
    
    % Calculate the derivative of h with respect to theta
    h_theta = DNQ_wrt_theta*e;
    
    % Build the matrix D_{k}(T+N)Q_{k}
    DTNQ = BuildDTQ(fw + z_fw , alpha(ite).*(gw + z_gw),k);
    
    % Calculate the paritial derivative of D_{k}(T+N)Q_{k} with respect to
    % alpha
    DTNQ_alpha = BuildDTQ(fw_wrt_alpha + Partial_zfw_wrt_alpha,...
        Partial_gw_wrt_alpha + Partial_zgw_wrt_alpha,k);
    
    % Calculate the paritial derivative of D_{k}(T+N)Q_{k} with respect to
    % theta
    
    DTNQ_theta = BuildDTQ(fw_wrt_theta + Partial_zfw_wrt_theta,...
        alpha(ite).*(gw_wrt_theta + Partial_zgw_wrt_theta) ...
        ,k);
    
    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    DY = BuildDY_Roots(m,n,k,idx_col,x_ls,alpha(ite),th(ite),ratio);
    
    % Calculate the matrix DP where P is the matrix such that c = P[f;g]
    DP  = BuildDP_Roots(m,n,alpha(ite),th(ite),idx_col,k,ratio);
    
    % Calculate the residual vector
    res_vec = (ct + h) - DTNQ*M*xk;
    
    % Create matrix E.
    E = eye(2*m+2*n-2*k+5);
    
    % Create the matrix C. This is made up of four submatrices, HZ, Hx,
    % H_alpha and H_theta.
    
    Hz      = GetH_z(idx_col,n,k,DY,DP);
    
    Hx      = DTNQ*M;
    
    H_alpha = DTNQ_alpha*M*xk - (ct_wrt_alpha+ h_alpha);
    
    H_theta = DTNQ_theta*M*xk - (ct_wrt_theta+ h_theta);
    
    C = [Hz,Hx,H_alpha,H_theta];  % the matrix C
    
    % Calculate the new right hand vector
    ek = ct + h;
    
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(res_vec) / norm(ek);
    
    % Update fnew - used in LSE Problem.
    p = -(yy-start_point);
    
    
    
end

if(SETTINGS.PLOT_GRAPHS)
    
    figure_name = sprintf('%s : SNTLN Residuals',mfilename);
    figure('name',figure_name)
    hold on
    title('Residuals over SNTLN iterations')
    plot(log10(condition),'-s')
    xlabel('iterations')
    ylabel('Residuals')
    hold off
    
end

if ite == SETTINGS.MAX_ITERATIONS_SNTLN
    fprintf('SNTLN_Roots : failed to converge after %i iterations \n\n',ite)
    fx_output = fx_n;
    gx_output = gx_new;
    X_output = x_ls;
    alpha_output = i_alpha;
    theta_output = theta;
    return
end

% Once iterations are complete, assign fx output, gx output, solution X
% output, alpha output and theta output.
fx_output = fx_n + zk(1:m+1);


for i = 0:1:length(z_fx)-2
    z_gx(i+1) = m*(z_fx(i+2) - z_fx(i+1))  ;
end

gx_output = gx_n + z_gx;

X_output  = xk;

alpha_output = alpha(ite);

theta_output = th(ite);

% Print the number of iterations

fprintf(' %s : Iterations required : %i \n', mfilename,ite);


figure_name = sprintf('%s : SNTLN Residuals',mfilename);
figure('name',figure_name)
title('Residual at each Iteration of SNTLN')
hold on
plot(1:1:length(condition),log10(condition));
xlabel('Iteration Number')

if(SETTINGS.PLOT_GRAPHS0
    
        figure_name = sprintf('%s : Condition',mfilename);
        figure('name',figure_name)
        title('Residual at each Iteration of SNTLN')
        hold on
        plot(1:1:length(condition),log10(condition));
        xlabel('Iteration Number')
        
        
        figure_name = sprintf('%s : Thetas',mfilename);
        figure('name',figure_name)
        title('Theta at each Iteration of SNTLN')
        hold on
        plot(1:1:length(th),log10(th));
        xlabel('Iteration Number')
        
        figure_name = sprintf('%s : Alphas',mfilename);
        figure('name',figure_name)
        title('Alpha at each iteration of SNTLN');
        hold on
        plot(1:1:length(alpha),log10(alpha));
end


end


function Hz = GetH_z(mincol,n,t,DY,DP)
if mincol<=(n-t+1)
    Hz = (DY-DP);
else
    Hz = DY-DP;
end
end







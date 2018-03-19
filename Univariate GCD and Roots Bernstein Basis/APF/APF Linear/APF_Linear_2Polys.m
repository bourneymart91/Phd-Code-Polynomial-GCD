function [fx_lr, gx_lr, dx_lr, ux_lr, vx_lr] = ...
    APF_Linear_2Polys(fx, gx, ux, vx, k)
% Refine Approximate Polynomial Factorisation by APF
%
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x) in the Bernstein basis
%
% ux : (Vector) Coefficients of the polynomial u(x), where v(x) is the
% quotient polynomial such that f(x)/u(x) = d(x)
% 
% vx : (Vector) Coefficients of the polynomial v(x), where v(x) is the
% quotient polynomial such that g(x)/v(x) = d(x)
%
% k : (Int) Calculated degree of d(x)
%
% % Outputs
%
% fx_lr : (Vector) Coefficients of the polynomial f(x), with added
% perturbations
%
% gx_lr : (Vector) Coefficients of the polynomial g(x), with added
% perturbations
%
% dx_lr : (Vector) Coefficients of the polynomial d(x), where d(x) is the
% greatest common divisor of f(x) and g(x)
%
% ux_lr : (Vector) Coefficients of the polynomial u(x), where v(x) is the
% quotient polynomial such that f(x)/u(x) = d(x), after perturbations have
% been added to f(x) and g(x)
%
% vx_lr : (Vector) Coefficients of the polynomial v(x), where v(x) is the
% quotient polynomial such that g(x)/v(x) = d(x), after perturbations have
% been added to f(x) and g(x)


% Global Variables

global SETTINGS

% Initialise iteration index
ite = 1;

% Initialise

% Get degree of polynomial f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Get number of coefficients in u(x), v(x), f(x) and g(x)
nCoefficients_ux = m - k + 1;
nCoefficients_vx = n - k + 1;
nCoefficients_fx = m + 1;
nCoefficients_gx = n + 1;



% Initialise S, the vector of thetas corresponding to coefficients of f(x),
% such that s_{k} = S * p_{k}
th_f = eye(m + 1);

% Initialise T - Matrix such that tk = T * qt
th_g = eye(n + 1);

% Initialise zk - Structured perturbations of u and v
zk = zeros(m+n-2*k+2,1);

% Get H^{-1} * C(u,v) * G
[HCG, ~, ~] = BuildHCG_2Polys(ux,vx,k);

%Build the RHS vector b = [fx; gx]
bk = [fx ; gx];

% Get d(x)
dx = SolveAx_b(HCG,bk);

% Get initial residual
res_vec = bk - ((HCG)*dx);

% Set some initial values
z_fx = zeros(nCoefficients_fx, 1);
z_gx = zeros(nCoefficients_gx, 1);

% Set initial values for the iterative process
z1_ux = zeros(nCoefficients_ux, 1);
z2_vx = zeros(nCoefficients_vx, 1);

% Construct the coefficient matrix in the equation that defines
% the constraint for the LSE problem.
HYk         = BuildHYQ_SNTLN(dx, m, n, 1);


% % Build the matrix C given by Hz Hp Hq

% Build the matrix H_{z}
H_z = HYk;

% Build the matrix H_{p}
H_p = -1.* th_f;

% Build the matrix H_{q}
H_q = -1.* th_g;

% Build the matrix C
C_temp      = ...
    [
    H_p,             zeros(m+1,n+1);...
    zeros(n+1,m+1),  H_q, ...
    ];

%
C = [H_z , C_temp];

%
E = eye(nCoefficients_fx + nCoefficients_gx + nCoefficients_ux + nCoefficients_vx);

%
ek = bk;

% Get the condition
condition(ite) = norm(res_vec)/norm(ek);

start_point = [...
        zk;...
        z_fx;...
        z_gx;...
    ];

%
yy = start_point;

%
f = -(yy - start_point);

% Start the iterative procedure for the solution of the LSE problem.
while condition(ite) > (SETTINGS.MAX_ERROR_APF) && ite < SETTINGS.MAX_ITERATIONS_APF
    
    % Use the QR decomposition to solve the LSE problem.
    % min |y-p| subject to Cy=q
    y = LSE(E, f, C, res_vec);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % % obtain the small changes.
    
    % Get change in z_{k} = [z_{u} z_{v}]
    delta_zk = y(1:nCoefficients_ux + nCoefficients_vx);
    
    % Get change in z_{f}(x)
    %     delta_zf_k = y(m+n-2*k+3:2*m+n-2*k+3);
    delta_zf_k = y(nCoefficients_ux + nCoefficients_vx + 1: nCoefficients_ux + nCoefficients_vx + nCoefficients_fx );
    
    % Get change in z_{g}(x)
    delta_zg_k = y(2*m+n-2*k+4 : 2*m+2*n-2*k+4);
        
    % Update variables zk, pk, qk,
    
    % Update z_{k}
    zk = zk + delta_zk;
    
    % Update z_{f}(x)
    z_fx = z_fx + delta_zf_k;
    
    % Update z_{g}(x)
    z_gx = z_gx + delta_zg_k;
    
       
    % Get perturbation vector, and separate in to perturbations of f,
    % z1 and perturbations of g, z2
    z1_ux = zk(1:m-k+1);
    z2_vx = zk(m-k+2:m+n-2*k+2);
    
    
    % % Build the components of the matrix C
    
    % Calculate H_z
    HYk = BuildHYQ_SNTLN(dx, m, n, 1);
    
    % Build H_z
    H_z = HYk;
    
    % Calculate H_p
    H_p = -1.*th_f;
    
    % Calculate H_q
    H_q = -1.*th_g;
    
    
    % % Build the matrix C
    C_temp = [  H_p,                zeros(m+1,n+1);
                zeros(n+1,m+1),     H_q, ];
    
    % Build Matrix C
    C = [H_z , C_temp];
    
    % Calculate Matrix H(C+E)G
    [HCEG,~,~] = BuildHCG_2Polys(ux + z1_ux, vx + z2_vx,k);
    
    % Calculate the new residual
    res_vec = [(fx + z_fx) ; (gx + z_gx)] - ((HCEG)*dx);
    
    % Calculate the new right hand vector.
    ek = [(fx + z_fx) ; (gx + z_gx)];
       
    % Update Condition scalar.
    condition(ite) = norm(res_vec)/norm(ek);
    
    % Update fnew
    f = -(yy - start_point);
    
    
end

% Display number of iterations
LineBreakLarge();
fprintf('Iterations over Nonlinear approximate polynomial factorisation : %i \n', ite);
LineBreakLarge();
SETTINGS.APF_REQ_ITE = ite;


% Update values of quotients u and v,
ux_lr = ux + z1_ux;
vx_lr = vx + z2_vx;

% Update value of common divisor dx
dx_lr = dx;

% Edit 20/07/2015
fx_lr = (fx + z_fx);
gx_lr = (gx + z_gx);


Plot_APF()
end





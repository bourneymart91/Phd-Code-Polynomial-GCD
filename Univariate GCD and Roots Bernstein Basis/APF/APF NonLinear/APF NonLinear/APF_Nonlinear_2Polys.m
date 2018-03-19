function [fx_lr, gx_lr, dx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = ...
    APF_Nonlinear_2Polys(fx, gx, ux, vx, i_alpha, i_theta, k)
% Refine Approximate Polynomial Factorisation by APF
%
%
% % Inputs
%
% fx : (Vector) Coefficients of the polynomial f(x)
% 
% gx : (Vector) Coefficients of the polynomial g(x) 
%
% ux : (Vector) Coefficients of the polynomial u(x), where u(x) is the
% quotient polynomial f(x)/u(x) = d(x)
%
% vx : (vector) Coefficients of the polynomial v(x), where v(x) is the
% quotient polynomial g(x)/v(x) = d(x)
%
% i_alpha : (Float) Initial value of alpha
%
% i_theta : (Float) Initial value of \theta
%
% k : (Int) Calculated degree of d(x)
%
% % Outputs
%
% fx_lr : (Vector) 
%
% gx_lr : (Vector)
%
% dx_lr : (Vector)
%
% ux_lr : (Vector) 
%
% vx_lr : (Vector)



% Global Variables

global SETTINGS

% Initialise iteration index
ite = 1;

alpha(1) = i_alpha;
theta(1) = i_theta;

% Get degree of polynomial f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Get number of coefficients in u(x)
nCoefficients_ux = m-k+1;
nCoefficients_vx = n-k+1;
nCoefficients_fx = m+1;
nCoefficients_gx = n+1;

% Initialise some useful vectors
vec_m    = (0:1:m)';
vec_n    = (0:1:n)';
vec_k    = (0:1:k)';

% Initialise zk - Structured perturbations of u(x) and v(x)
zk = zeros(nCoefficients_ux + nCoefficients_vx,1);

% Convert f and g to modified bernstein basis, excluding binomial
% coefficient
fw = GetWithThetas(fx, theta);
gw = GetWithThetas(gx, theta);

% Convert u and v to modified bernstein basis, excluding binomial
% coefficient
uw = GetWithThetas(ux, theta);
vw = GetWithThetas(vx, theta);

% Initialise S = th_f, the vector of thetas corresponding to coefficients
% of f(x), such that s_{k} = S * p_{k}
th_f = (diag(theta(ite).^vec_m));

% Initialise T - Matrix such that tk = T * qt
th_g = (diag(theta(ite).^vec_n));

% Initialise zk - Structured perturbations of u and v
zk = zeros(m+n-2*k+2,1);

% Get partial derivatives of f(\omega) and g(\omega) with respect to theta
fw_wrt_theta = Differentiate_wrt_theta(fw, theta(ite));
gw_wrt_theta = Differentiate_wrt_theta(gw, theta(ite));

% Get partial derivatives of uw and vw with respect to theta.
Partial_uw_wrt_theta = Differentiate_wrt_theta(uw, theta(ite));
Partial_vw_wrt_theta = Differentiate_wrt_theta(vw, theta(ite));

% Get H^{-1} * C(u,v) * G
[HCG, H1C1G, H2C2G] = BuildHCG_2Polys(uw, vw, k);


% Build HCG with respect to theta
[~,H1C1G_wrt_theta, H2C2G_wrt_theta] = ...
    BuildHCG_2Polys(Partial_uw_wrt_theta, Partial_vw_wrt_theta, k);

%Build the RHS vector b = [fw,alpha.*gw]
bk = [fw ; alpha(ite).*gw];

dw = SolveAx_b(HCG, bk);
dx = GetWithoutThetas(dw, theta(1));

% Get partial derivative of dw with respect to theta
dw_wrt_theta = Differentiate_wrt_theta(dw, theta(ite));

% Get initial residual
res_vec = bk - ((HCG)*dw);

% Set some initial values
residual(ite)   = norm(res_vec);
z_fx = zeros(nCoefficients_fx, 1);
z_gx = zeros(nCoefficients_gx, 1);

% Values of LHS Perturbations
% Obtain structured perturbations sw of fw, and tw of gw
z_fw = GetWithThetas(z_fx, theta);
z_gw = GetWithThetas(z_gx, theta);

% Set initial values for the iterative process
z1_ux = zeros(length(uw), 1);
z2_vx = zeros(length(vw), 1);

% Construct the coefficient matrix in the equation that defines
% the constraint for the LSE problem.
HYk         = BuildHYQ_SNTLN(dx, m, n, theta(ite));


% % Build the matrix C given by Hz Hp Hq Halpha Htheta1 Htheta2
H_z         = HYk;

H_p         = (-1)*th_f;

H_q         = -(alpha(ite))*th_g;

H_alpha     = -(gw + z_gw);

H_theta1    = -fw_wrt_theta+...
    (H1C1G_wrt_theta * dw)+...
    (H1C1G * dw_wrt_theta);

H_theta2    = (-(alpha(ite))*gw_wrt_theta)+...
    (H2C2G_wrt_theta * dw)+...
    (H2C2G * dw_wrt_theta);

C_temp      = ...
    [
    H_p,             zeros(m+1,n+1), zeros(m+1,1), H_theta1;...
    zeros(n+1,m+1),  H_q,            H_alpha,       H_theta2...
    ];

C       = [H_z , C_temp];

E       = eye(2*m+2*n-2*k+6);

ek = bk;

% Get the condition
condition(ite) = norm(res_vec)/ norm(ek);

start_point = [...
    zk;...
    z_fx;...
    z_gx;...
    alpha;...
    theta];

yy = start_point;

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
    delta_zg_k = y(2*m+n-2*k+4:2*m+2*n-2*k+4);
    
    % Get change in \alpha
    delta_alpha     = y(end-1);
    
    % Get change in \theta
    delta_theta     = y(end);
    
    % Update variables zk, pk, qk, beta, theta
    
    % Update z_{k}
    zk = zk + delta_zk;
    
    % Update z_{f}(x)
    z_fx = z_fx  + delta_zf_k;
    
    % Update z_{g}(x)
    z_gx = z_gx  + delta_zg_k;
    
    % Update \alpha
    alpha(ite) = alpha(ite-1) + delta_alpha;
    
    % Update \theta
    theta(ite) = theta(ite-1) + delta_theta;
    
    % %
    % Update the polynomial f(\omega)
    fw = GetWithThetas(fx, theta(ite));
    
    % Update the polynomial g(\omega)
    gw = GetWithThetas(gx, theta(ite));
    
    % Update matrices S = th_f and T = th_g
    th_f = diag(theta(ite).^vec_m);
    th_g = diag(theta(ite).^vec_n);
    
    % Update the polynomial d(\omega)
    dw = GetWithThetas(dx, theta(ite));
    
    % Update the polynomial u(\omega)
    uw = GetWithThetas(ux, theta(ite));
    
    % Update the polynomial v(\omega)
    vw = GetWithThetas(vx, theta(ite));
    
    % Obtain partial derivatives of u(\omega) and v(\omega) with respect to
    % \theta
    uw_wrt_theta = Differentiate_wrt_theta(uw, theta(ite));
    vw_wrt_theta = Differentiate_wrt_theta(vw, theta(ite));
    
    % Build Matrices H_{1}C_{1}(u)G and H_{2}C_{2}(v)G
    [~, H1C1G, H2C2G] = BuildHCG_2Polys(uw, vw, k);
    
    % Build Matrices H_{1}C_{1}(u)G and H_{2}C_{2}(v)G with respect to
    % theta
    [~, H1C1G_wrt_theta, H2C2G_wrt_theta] = BuildHCG_2Polys(uw_wrt_theta, vw_wrt_theta,k);
    
    % Get perturbation vector, and separate in to perturbations of f,
    % z1 and perturbations of g, z2
    z1_ux = zk(1:m-k+1);
    z2_vx = zk(m-k+2:m+n-2*k+2);
    
    % Obtain z1 and z2 in the modified bernstein basis.
    z_uw = GetWithThetas(z1_ux, theta(ite));
    z_vw = GetWithThetas(z2_vx, theta(ite));
    
    % obtain partial derivatives of z1 and z2 with respect to theta
    z2w_wrt_theta = Differentiate_wrt_theta(z_vw, theta(ite));
    z1w_wrt_theta = Differentiate_wrt_theta(z_uw, theta(ite));
    
    % Build Matrices H_{1}E_{1}(z1)G and H_{2}E_{2}(z2)G
    [~, H1E1G, H2E2G] = BuildHCG_2Polys(z_uw, z_vw, k);
    
    % Calculate Partial derivatives of Matrices H_{1}E_{1}(z1)G and
    % H_{2}E_{2}(z2)G with respect to theta
    [~, H1E1G_wrt_theta, H2E2G_wrt_theta] = BuildHCG_2Polys(z1w_wrt_theta, z2w_wrt_theta, k);
    
    % Obtain structured perturbations sw of fw, and tw of gw
    z_fw = GetWithThetas(z_fx, theta(ite));
    z_gw = GetWithThetas(z_gx, theta(ite));
    
    % Calculate partial derivatives of sw and tw with respect to theta
    s_wrt_theta = Differentiate_wrt_theta(z_fw, theta(ite));
    t_wrt_theta = Differentiate_wrt_theta(z_gw, theta(ite));
    
    % Calculate partial derivatives of fw and gw with respect to theta
    fw_wrt_theta = Differentiate_wrt_theta(fw, theta(ite));
    gw_wrt_theta = Differentiate_wrt_theta(gx, theta(ite));
    
    % Calculate partial derivative of dw with respect to theta
    dw_wrt_theta = Differentiate_wrt_theta(dw, theta(ite));
    
    % Build Matrix C
    
    % Calculate H_z
    HYk = BuildHYQ_SNTLN(dx, m, n, theta(ite));
    
    % Build H_z
    H_z = HYk;
    
    % Calculate H_p
    H_p = (-1)* th_f;
    
    % Calculate H_q
    H_q = -(alpha(ite))* th_g;
    
    % Calculate H_beta
    H_alpha = -(gw + z_gw);
    
    % Calculate H_theta1
    H_theta1 = -(fw_wrt_theta + s_wrt_theta)+...
        ((H1C1G_wrt_theta)* dw)+...
        ((H1E1G_wrt_theta)* dw)+...
        ((H1C1G)* dw_wrt_theta)+...
        ((H1E1G)* dw_wrt_theta);
    
    % Calculate H_theta2
    H_theta2 = -(alpha(ite))*(gw_wrt_theta + ...
        t_wrt_theta)+...
        ((H2C2G_wrt_theta)*dw)+...
        ((H2E2G_wrt_theta)*dw)+...
        ((H2C2G)*dw_wrt_theta)+...
        ((H2E2G)*dw_wrt_theta);
    
    C_temp = [H_p, zeros(m+1,n+1), zeros(m+1,1), H_theta1;
        zeros(n+1,m+1), H_q, H_alpha      , H_theta2];
    
    % Build Matrix C
    C = [H_z , C_temp];
    
    % Calculate Matrix H(C+E)G
    [HCEG,~,~] = BuildHCG_2Polys(uw+z_uw,vw+z_vw,k);
    
    % Calculate the new residual
    res_vec = [fw + z_fw ;(alpha(ite)*(gw + z_gw))]-((HCEG)*dw);
    
    % Calculate the new right hand vector.
    ek = ...
        [
        fw + z_fw ;...
        (alpha(ite)*(gw + z_gw))...
        ];
    
    
    % update vector of residual
    residual(ite) = norm(res_vec);
    
    % Update Condition scalar.
    condition(ite) = norm(res_vec)/norm(ek);
    
    % Update fnew
    f = -(yy - start_point);
    
    % Edit 01/06/2015 16:30:00
    
    % Calculate the termination criterion in the modified
    % Bernstein basis..
    [res_uw(ite), res_vw(ite)] = Term_Criterion_APF(fw,gw,z_fw,z_gw,uw,vw,...
        dw,k,alpha(ite));
    
    % Repeat this calculation for the Bernstein basis. Transform the
    % variables from the modified Bernstein basis to the Bernstein
    % basis.
    
    fx_p = GetWithoutThetas(fw, theta(ite));
    sx_p = GetWithoutThetas(z_fw, theta(ite));
    gx_p = GetWithoutThetas(gw, theta(ite));
    tx_p = GetWithoutThetas(z_gw, theta(ite));
    ukx = ux + z1_ux;
    vkx = vx + z2_vx;
    dkx = dw./(theta(ite).^vec_k);
    
    [res_ux(ite),res_vx(ite)] = Term_Criterion_APF(fx_p,gx_p,sx_p,...
        tx_p, ukx, vkx, dkx, k, 1.0);
    
    
end

% Plot graphs
Plot_APF()

% Display number of iterations
LineBreakLarge();
fprintf('Iterations over approximate polynomial factorisation : %i \n', ite+1);
LineBreakLarge();
SETTINGS.APF_REQ_ITE = ite;


if(SETTINGS.PLOT_GRAPHS)
    
    
    % Plot the normalised residuals res_ux, res_vx, res_uw and res_vw.
    if(ite > 1)
        %plotgraphs3(res_ux, res_vx, res_uw, res_vw);
        % Write out the number of iterations required and plot the values of
        % alpha, theta and the residual.
        %plotgraphs4(alpha, theta, residual);
    end
    
    
end



% Update values of quotients u and v,
ux_lr = ux + z1_ux;
vx_lr = vx + z2_vx;

% Update value of theta
theta_lr = theta(ite);
alpha_lr = alpha(ite);

% Update value of common divisor dx
dx_lr = GetWithoutThetas(dw, theta(ite));

% Edit 20/07/2015
fx_lr = GetWithoutThetas((fw + z_fw), theta(ite));
gx_lr = GetWithoutThetas((gw + z_gw), theta(ite));




end





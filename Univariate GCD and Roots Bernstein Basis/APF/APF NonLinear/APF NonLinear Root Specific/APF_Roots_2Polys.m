function [PostAPF_fx, PostAPF_gx, PostAPF_dx, PostAPF_uk, PostAPF_vk, PostAPF_theta] = ...
    APF_Roots_2Polys(fx,ux,vx,i_theta,dx,t)
% APF_Roots(fx,ux,vx,initial_theta,dx,t)
%
% Refine Approximate Polynomial Factorisation by newton raphson method iteration
%
%
% % Inputs
%
% fx : Coefficients of the polynomial f(x) in the Bernstein basis
%
% gx : Coefficients of the polynomial g(x) in the Bernstein basis
%
% ux : Coefficients of the polynomial u(x) in the Bernstein basis
%
% vx : Coefficients of the polynomial v(x) in the Bernstein basis
%
% alpha - Optimal value of \alpha
%
% i_theta : Optimal value of \theta
%
% dx : initial approximation of GCD(fx,gx)
%
% t : Calculated degree of GCD
%
% % Outputs
%
%

global SETTINGS


% Initialise iteration index
ite = 1;

% Set initial value of \theta.
th(ite) = i_theta;

% Set the initial value of \alpha.
%alpha(ite) = initial_alpha;

% Get degree of polynomial f(x)
m = GetDegree(fx);

% Initialise some useful vectors
vecm  = (0:1:m)';

% Get coefficients of the polynomial f(\omega)
fw_n = GetWithThetas(fx,th(ite));

% Get coefficients of the polynomial g(\omega)
uw = GetWithThetas(ux,th(ite));

% Get coefficients of the polynomail d(\omega)
dw = GetWithThetas(dx,th(ite));

% Initialise S, such that sk = S * p_{t}, where p_{t} is a vector of
% perturbations applied to u
S = (diag(th(ite).^vecm));

% Initialise z_{u}(x), the perturbations of u(x)
v_zu = zeros(m-t+1,1);

% Get partial derivatives of f(\omega) and g(\omega) with respect to \theta
fw_wrt_theta = Differentiate_wrt_theta(fw_n,th(ite));

% Get partial derivative of d(\omega) with respect to \theta
dw_wrt_theta = Differentiate_wrt_theta(dw,th(ite));

% Get partial derivatives of u(\omega) and v(\omega) with respect to \theta
uw_wrt_theta = Differentiate_wrt_theta(uw,th(ite));

% Build H_{t}C(f,g)G_{t}
H1C1G = BuildH1C1G(uw,t);

% Build HCG with respect to theta
H1C1G_wrt_theta = BuildH1C1G(uw_wrt_theta,t);

%Build the RHS vector b = [fw,alpha.*gw]
bk = fw_n;

% Get initial residual vector
res_vec = bk - ((H1C1G)*dw);

% %
% Set some initial values

% The initial residual vector
vResidual(ite) = norm(res_vec);

% Initialise vector z_{f}(x), perturbations to polynomial f(x)
v_zf               = zeros(m+1,1);

% Initialise vector z_{u}(x), perturbations of polynomial u(x)
v_zu = zeros(length(uw),1);

% %
% Construct the coefficient matrix in the equation that defines
% the constraint for the LSE problem.

HYk        = BuildHYQ_Roots(dx,m,th(ite));

H_z         = HYk;

H_p         = (-1)*S;

H_theta1    = -fw_wrt_theta+...
    (H1C1G_wrt_theta * dw)+...
    (H1C1G * dw_wrt_theta);

C_temp      = [ H_p,            H_theta1];

C       = [H_z , C_temp];

E       = eye(2*m-t+3);

fnew    = zeros(2*m-t+3,1);

ek = bk;

condition = norm(res_vec)/norm(ek);

startpoint = [v_zu;v_zf;th];

yy = startpoint;

% Start the iterative procedure for the solution of the LSE problem.

while condition > (SETTINGS.MAX_ERROR_APF) && ite < SETTINGS.MAX_ITERATIONS_APF
    
    % Use the QR decomposition to solve the LSE problem.
    % min |y-p| subject to Cy=q
    y = LSE(E,fnew,C,res_vec);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % %
    % Obtain the small changes
    delta_zu    = y(1:m-t+1);
    delta_zf    = y(m-t+2:2*m-t+2);
    delta_theta = y(end);
    
    % %
    % Update
    
    % Update the vector z_{u}(x)
    v_zu = v_zu + delta_zu;
    
    % Update the vector of z_{f}(x)
    v_zf = v_zf  + delta_zf;
    
    th(ite)  = th(ite-1) + delta_theta;
    
    % Update the iterative value of f and g
    fw = GetWithThetas(fx,th(ite));
    
    % Update matrices S and T
    S = diag(th(ite).^vecm);
    
    % Update the iterative value of GCD dw
    dw = GetWithThetas(dx,th(ite));
    
    % Update the iterative value of quotients uw and vw
    uw = GetWithThetas(ux,th(ite));
    
    % Obtain partial derivatives of uw and vw with respect to theta
    uw_wrt_theta = Differentiate_wrt_theta(uw,th(ite));
    
    % Build Matrices H_{1}C_{1}(u)G and H_{2}C_{2}(v)G
    H1C1G = BuildH1C1G(uw,t);
    
    % Build Matrices H_{1}C_{1}(u)G and H_{2}C_{2}(v)G with respect to
    % theta
    H1C1G_wrt_theta = BuildH1C1G(uw_wrt_theta,t);
    
    % Get perturbations of polynomial f
    v_zu = v_zu(1:m-t+1);
    
    % Obtain z1 and z2 in the modified bernstein basis.
    z1w = GetWithThetas(v_zu,th(ite));
    
    % obtain partial derivatives of z1 and z2 with respect to theta
    z1w_wrt_theta = Differentiate_wrt_theta(z1w,th(ite));
    
    % Build the matrix H1*E1*G
    H1E1G = BuildH1C1G(z1w,t);
    
    % Calculate Partial derivatives of Matrices H_{1}E_{1}(z1)G and
    % H_{2}E_{2}(z2)G with respect to theta
    H1E1G_wrt_theta = BuildH1C1G(z1w_wrt_theta,t);
    
    % Obtain structured perturbations sw of fw, and tw of gw
    sw   = GetWithThetas(v_zf,th(ite));
    
    % Calculate partial derivatives of sw and tw with respect to theta
    s_wrt_theta = Differentiate_wrt_theta(sw,th(ite));
    
    % Calculate partial derivatives of fw and gw with respect to theta
    fw_wrt_theta    = Differentiate_wrt_theta(fw,th(ite));
    
    % Calculate partial derivative of dw with respect to theta
    dw_wrt_theta    = Differentiate_wrt_theta(dw,th(ite));
    
    % Calculate H_z
    
    HYk         = BuildHYQ_Roots(dx,m,th(ite));
    
    H_z = HYk;
    
    % Calculate H_p
    H_p = (-1)*S;
    
    % Calculate H_theta1
    H_theta1 = -(fw_wrt_theta + s_wrt_theta)+...
        ((H1C1G_wrt_theta)*dw)+...
        ((H1E1G_wrt_theta)*dw)+...
        ((H1C1G)*dw_wrt_theta)+...
        ((H1E1G)*dw_wrt_theta);
    
    C_temp = [  H_p, H_theta1];
    
    % Build Matrix C
    C = [H_z , C_temp];
    
    % Build the matrix HCEG
    HCEG = BuildH1C1G(uw+z1w,t);
    
    % Calculate the new residual
    res_vec = fw + sw - ((HCEG)*dw);
    
    % Calculate the new right hand vector.
    ek = fw + sw;
    
    % update vector of residual
    vResidual(ite) = norm(res_vec);
    
    % Update Condition scalar.
    condition = norm(res_vec)/norm(ek);
    
    % Update fnew
    fnew = -(yy - startpoint);
    
    
end

if ite == SETTINGS.MAX_ITERATIONS_APF
    fprintf('APF failed to converge after %i iterations\n',ite)
    PostAPF_dx = dx;
    PostAPF_uk = ux;
    PostAPF_vk = vx;
    PostAPF_theta = initial_theta;
    
    % Edit 20/07/2015
    PostAPF_fx = fx;
    PostAPF_gx = gx;
    return
end


% Display number of iterations
LineBreakLarge()
fprintf('Iterations over approximate polynomial factorisation : %i \n', ite+1);
LineBreakLarge()

% Update values of quotients u(x) and v(x)
PostAPF_uk = ux + v_zu;

PostAPF_vk = vx;

% Update value of \theta
PostAPF_theta = th(ite);

% Update value of common divisor dx
PostAPF_dx = GetWithoutThetas(dw,th(ite));

% update value of the input polynomial f
PostAPF_fx = GetWithoutThetas(fw + sw,th(ite));

% update values of the input polynomial g, which is now calculated by
% derivative constraint.
PostAPF_gx = zeros(m,1);
% for i = 0:1:m-1
%    PostAPF_gx(i+1) = m.*(PostAPF_fx(i+2) - PostAPF_fx(i+1));
% end


if(SETTINGS.PLOT_GRAPHS)
    figure_name = sprintf([mfilename 'Plot Residuals']);
    figure('name',figure_name)
    title('plotting residual')
    hold on
    plot(1:1:length(vResidual),(vResidual));
    hold off
end

end










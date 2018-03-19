function [fx_lr, gx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = STLN(fx,gx,k)
% Perform Structured Total Least Norm to obtain a low rank approximation
% of the t-th Sylvester matrix. Note this is a linear problem, any
% alpha and theta values are already included in f(x) and g(x).
%
%
% % Inputs.
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% gx : (Vector) Coefficients of polynomial g(x)
%
% t : (Int) Degree of GCD(f(x),g(x))
%
% % Outputs.
%
% fx_lr : (Vector) Coefficients of f(x) after refinement f(x) + \delta
%
% gx_lr : (Vector) Coefficients of g(x) after refinement g(x) + \delta
%
% 
%
global SETTINGS

% Get degree of polynomial f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Initialise the vector of perturbations zf(x)
zf = zeros(m+1,1);

% Initialise the vector of perturbations zg(x)
zg = zeros(n+1,1);

z = [zf ; zg];

% Build the k'th subresultant
T1_fx = BuildT1(fx, n-k);
T1_gx = BuildT1(gx, m-k);
Sk_fg = [T1_fx T1_gx];

% Build the matrix E_{t}(z)
T1_zf = BuildT1(zf, n-k);
T2_zg = BuildT1(zg, m-k);
Sk_zfzg = [T1_zf T2_zg];

% Get the index of the optimal colummn for removal
[~,idx_col] = GetMinDistance(Sk_fg);

% Get A_{t} the LHS matrix, equivalent to S_{t} with the optimal column
% removed.
At_fg = Sk_fg;
At_fg(:,idx_col) = [];

% Get c_{t} the removed column of S_{t} to form A_{t}.
ck = Sk_fg(:,idx_col);

% Get E_{t}, the matrix of structured perturbations corresponding to A_{t}.
Ak_zfzg = Sk_zfzg;
Ak_zfzg(:,idx_col) = [];

% Get h_{t}, the vector of structured perturbations corresponding to c_{t}
hk = Sk_zfzg(:,idx_col);

% Build Pt
Pk = BuildP_STLN(m, n, k, idx_col);

% Get the solution vector x of A_{t}x = c_{t}.
xk = SolveAx_b(At_fg,ck);

% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
res_vec = (ck + hk) - At_fg*xk;

% Build the matrix Y_{t}
x = [xk(1:idx_col-1) ; 0 ; xk(idx_col:end)];

Yk = BuildY_STLN(x, m, n, k);

% Build the matrix C for LSE Problem
H_z = Yk - Pk;
H_x = At_fg + Ak_zfzg;

C = [H_z H_x];

% Build the matrix E for LSE Problem
E = eye(2*m+2*n-2*k+3);


% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
        z;...
        xk;
    ];

% Set yy to be the vector which stores all cummulative perturbations.
yy = start_point;

% Set the initial value of vector p to be zero
f = -(yy - start_point);


% Initialise the iteration counter
ite = 1;

% Set the termination criterion
condition(ite) = norm(res_vec)./ norm(ck);



while condition(ite) >  SETTINGS.MAX_ERROR_SNTLN &&  ite < SETTINGS.MAX_ITE_SNTLN
    
    % Increment interation counter
    ite = ite + 1;
    
    % Get small petrubations by LSE
    y_lse = LSE(E, f, C, res_vec);
    
    % Increment cummulative peturbations
    yy = yy + y_lse;
    
    % obtain the small changes to z and x
    delta_zk        = y_lse(1:m+n+2,1);
    
    % Update z and x
    z = z + delta_zk;
    
    % Split z into z_f and z_g
    zf = z(1:m+1);
    zg = z(m+2:end);
    
    % Build the matrix B_{t} = [E1(zf) E2(zg)]
    Sk_zfzg = BuildT(zf,zg,k);
    
    % Get the matrix E_{t} with optimal column removed
    Ak_zfzg = Sk_zfzg;
    Ak_zfzg(:,idx_col) = [];
    
    % Get the column vector h_{t}, the optimal column removed from B_{t},
    % and equivalent to c_{t} removed from S_{t}
    hk = Sk_zfzg(:,idx_col);
    
    
    delta_xk        = y_lse((m+n+3):(2*m+2*n-2*k+3),1);
    xk = xk + delta_xk;
    x = [xk(1:idx_col-1) ; 0 ; xk(idx_col:end)];
    
    % Build the matrix Y_{t} where Y_{t}(x)*z = E_{t}(z) * x
    Yk = BuildY_STLN(x,m,n,k);
    
    % Get the residual vector
    res_vec = (ck+hk) - ((At_fg+Ak_zfzg)*xk);
    
    % Update the matrix C
    H_z = Yk - Pk;
    H_x = At_fg + Ak_zfzg;
    C = [H_z H_x];
    
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    % Update the termination criterion
    condition(ite) = norm(res_vec)./norm(ck + hk) ;
    
end

Plot_SNTLN();

LineBreakLarge()
fprintf([mfilename ' : ' sprintf('Required number of iterations : %i \n',ite)])
LineBreakLarge()
SETTINGS.LOW_RANK_APPROX_REQ_ITE = ite;

% If the final condition is less than the original, output the new values,
% otherwise output the old values for f(x) and g(x).


fx_lr = fx + zf;
gx_lr = gx + zg;


% Get u(x) and v(x)

% Get the vector x
x = [xk(1:idx_col-1) ; -1 ; xk(idx_col:end)];

% Get the number of coefficients in v(x)
nCoefficients_vx = n-k+1;

% Get coefficients of v(x)
vx_lr = x(1:nCoefficients_vx);

% Get coefficients of u(x)
ux_lr = -1 .* x(nCoefficients_vx + 1 : end);





end



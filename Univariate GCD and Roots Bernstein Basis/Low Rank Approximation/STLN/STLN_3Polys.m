function [fx_lr, gx_lr, hx_lr, ux_lr, vx_lr, wx_lr] = STLN_3Polys(fx, gx, hx, k, idx_OptCol)
% Perform STLN with no refinement of alpha or theta to compute the low rank
% approximation of the Sylvester subresultant matrix S_{k}(f,g)
%
% Inputs.
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% gx : (Vector) coefficients of polynomial g(x)
%
% hx : (Vector) Coefficients of the polynomial h(x)
%
% k : (Int) Degree of GCD d(x)
%
% idx_col : (Int) Index of optimal column to be removed from S_{t}(f,g)
%
% Outputs.
%
% fx_lr : (Vector) Coefficients of f(x) after addition of structured perturbations
%
% gx_lr : (Vector) Coefficients of g(x) after addition of strucutred perturbations
%
% hx_lr : (Vector) Coefficients of the polynomial h(x)
%
% ux_lr : (Vector) Coefficients of u(x)
%
% vx_lr : (Vector) Coefficients of v(x)
%
% wx_lr : (Vector) Coefficients of w(x)

global SETTINGS

% Get the degree of polynomial f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

% Initialise the vector of perturbations zf(x)
z_fx = zeros(m + 1, 1);
z_gx = zeros(n + 1, 1);
z_hx = zeros(o + 1, 1);

% Initialise the vector of perturbations z.
z = [z_fx ; z_gx; z_hx];

% Build the t'th subresultant S_{t}(f,g)
D = BuildD_3Polys(m, n - k, o - k);
T = BuildT_3Polys(fx, gx, hx, k);
Q = BuildQ_3Polys(m, n, o, k);
DTQ = D * T * Q;

% Build the matrix E_{t}(z)
B = BuildT_3Polys(z_fx, z_gx, z_hx, k);
DBQ = D * B * Q;


% Get A_{t} the LHS matrix, equivalent to S_{t} with the optimal column
% removed.
Ak = DTQ;
Ak(:, idx_OptCol) = [];

% Get c_{t} the removed column of S_{t} to form A_{t}.
ck = DTQ(:, idx_OptCol);

% Get E_{t}, the matrix of strucured perturbations corresponding to A_{t}.
Bk = DBQ;
Bk(:, idx_OptCol) = [];

% Get h_{t}, the vector of strucutred perturbations corresponding to c_{t}
hk = DBQ(:, idx_OptCol);

% Build DP
DPQ = BuildDPG_STLN_3Polys(m, n, o, k, idx_OptCol);

% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
xk = SolveAx_b(Ak + Bk, ck + hk);

% Get residual vector
q = (ck + hk) - (Ak + Bk)*xk;

% Get the vector x with a zero included in the x_ls solution.
x = [xk(1 : idx_OptCol - 1) ; 0 ; xk(idx_OptCol : end)];

% Build the matrix Y_{t}
DYQ = BuildDYQ_STLN_3Polys(x, m, n, o, k);

% Build the Matrx C for LSE problem
C_z = DYQ - DPQ;
C_x = Ak + Bk;

C = [C_z C_x];

H1 =  eye(m + 1) .* (n - k + 1) .* (o - k + 1);
H2 =  eye(n + 1) .* (m - k + 1);
H3 =  eye(o + 1) .* (m - k + 1);

H = blkdiag(H1, H2, H3);
zeroblock = zeros(m + n + o + 3, m + n + o - (3*k) + 2);
E = [H zeroblock];



% Define the starting vector for the iterations for the LSE problem.
start_point = ...
    [...
        z;...
        xk;...
    ];


% Initialise yy the vector of accumulated perturbations.
yy = start_point;

% Set the initial value of vector p to be zero
p = (E * (start_point - yy));


% Initialise the iteration counter
ite = 1;

% Set the termination criterion.
vCondition(ite) = norm(q)./norm(ck);


while vCondition(ite) >  SETTINGS.MAX_ERROR_SNTLN &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN
    
    % increment iteration number
    ite = ite + 1;
    
    % Get small changes in vector y
    %y_lse = LSE(H, p, C, q);
    
    % usage: x = lse(A,b,C,d)
    % Minimizes norm(A*x - b),
    % subject to C*x = d
    y_lse = LSE(E, p, C, q);
   
    %y_lse = LSE(E,p,C,q);
    
    
    % add small changes to cummulative changes
    yy = yy + y_lse;
    
    % obtain the small changes in z and x
    delta_zk        = y_lse(1 : m + n + o + 3, 1);
    delta_xk        = y_lse((m + n + o + 4) : end, 1);
    
    % Update z and x
    z = z + delta_zk;
    xk = xk + delta_xk;
    
    
    % Split z into z_f and z_g
    z_fx = z(1 : m + 1);
    z_gx = z(m + 2  : m + n + 2);
    z_hx = z(m + n + 3 : end);
    
    
    % Build the matrix E = D^{-1} * B(zf,zg) * Q
    DBQ = D * BuildT_3Polys(z_fx, z_gx, z_hx, k) * Q;
    
    % Build the matrix Bt = DBQ with opt column removed.
    Bk = DBQ;
    Bk(:, idx_OptCol) = [];
    
    % Get h_{t}, the optimal column removed from BDQ
    hk = DBQ(:, idx_OptCol);
    
   
    % Get the vector x
    x = [xk(1 : idx_OptCol - 1) ; 0 ; xk(idx_OptCol : end)];
    
    % Separate the component parts of x into x_v and x_u, where x_v is an
    % approximation of v(x) and x_u is an approximation u(x).
    DYQ = BuildDYQ_STLN_3Polys(x, m, n, o, k);
          
    % Get updated residual vector
    q = (ck + hk) - ((Ak + Bk)*xk);
    
    % Update the matrix C
    C_z = DYQ - DPQ;
    C_x = Ak + Bk;
    C = [C_z C_x];
    
    % Update p - used in LSE Problem.
    p = E * (start_point - yy);
    
    % Update the termination criterion.
    vCondition(ite) = norm(q) ./ norm(ck + hk) ;
    
    
    
end

cond(Ak)
cond(Ak + Bk)


if SETTINGS.PLOT_GRAPHS_LOW_RANK_APPROXIMATION
    PlotResiduals(vCondition)
end





% Print the number of iterations required
LineBreakLarge()
fprintf([mfilename ' : ' sprintf('Required number of iterations : %i \n',ite)]);
LineBreakLarge()
SETTINGS.LOW_RANK_APPROX_REQ_ITE = ite;

% %
% Get polynomials for output.

% Get polynomial f(x) + \delta f(x)
fx_lr = fx + z_fx;
gx_lr = gx + z_gx;
hx_lr = hx + z_hx;

% % Get polynomials u(x) and v(x)

% Get u(x) and v(x) from x_ls
x_bar = [xk(1 : idx_OptCol - 1) ; -1 ; xk(idx_OptCol : end)];

% Get the number of coefficients in v(x)
nCoefficients_vx = n - k + 1;
nCoefficients_wx = o - k + 1;

% Get the polynomial v(x) from the vector x
vx_lr = x_bar(1 : nCoefficients_vx);
wx_lr = x_bar(nCoefficients_vx + 1 : nCoefficients_vx + nCoefficients_wx);

% Get the polynomial u(x) from the vector x
ux_lr = -1.* x_bar(nCoefficients_vx + nCoefficients_wx + 1 : end);


end



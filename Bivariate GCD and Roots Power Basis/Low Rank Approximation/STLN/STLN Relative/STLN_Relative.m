function [fxy_lr,gxy_lr,uxy_lr,vxy_lr] = STLN_Relative(fxy, gxy, k1, k2, idx_col)
% STLN_Relative(fxy, gxy, k1, k2, idx_col)
%
% Given coefficients f(x,y) and g(x,y) find the low rank approximation of
% the Syvlester subresultant S_{t_{1},t_{2}}.
%
%
% % Inputs.
%
% [fxy, gxy] : Coefficients of polynomial f(x,y) and g(x,y)
%
% k1 : Degree of d(x,y) with respect to x
%
% k2 : Degree of d(x,y) with respect to y
%
% idx_col : Index of optimal column for removal from S_{k_{1},k_{2}}(f,g)
%
% % Outputs
%
% fxy_lr : Coefficients of f(x,y) with added perturbations
%
% gxy_lr : Coefficients of g(x,y) with added perturbations
%
% uxy_lr : Coefficients of u(x,y)
%
% vxy_lr : Coefficients of v(x,y)

global SETTINGS

% Get degree of polynomials f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degree of polynomial g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Get the number of coefficients in the polynomial f(x,y)
nCoeffs_fxy = (m1+1) * (m2+1);

% Get the number of coefficients in the polynomial g(x,y)
nCoeffs_gxy = (n1+1) * (n2+1);

% Get the number of coefficients in v(x,y)
nCoeffs_vxy = (n1-k1+1) * (n2-k2+1);

% Get the number of coefficients in u(x,y)
nCoeffs_uxy = (m1-k1+1) * (m2-k2+1);

% Build the Sylvester Matrix S_{k1,k2}(f,g)
Sk_fg = BuildT_Relative_Bivariate_2Polys(fxy,gxy,k1,k2);


% Remove optimal column
Ak_fg = Sk_fg;
Ak_fg(:,idx_col) = [];
ck = Sk_fg(:,idx_col);

% Build the matrix A_{k1,k2}(zf,zg)
Ak_zfzg = zeros(size(Ak_fg));

% Build the vector removed from S_{k}(zf,zg)
hk = zeros(size(ck));

% Initialise vector of perturbations of f and g
z = zeros(nCoeffs_fxy + nCoeffs_gxy,1);

% Get the vector of coefficients zf(x,y)
v_zf = z(1:nCoeffs_fxy );

% Get the vector of coefficeints zg
v_zg = z(nCoeffs_fxy+1:end);

% Get zf as a matrix
% EDIT 10/03/2016 - vZ_fxy has zeros removed, so include the zeros to form
% a matrix m1+1 * m2+1
zf_xy = GetAsMatrix(v_zf, m1, m2);

% Get zg as a matrix
zg_xy = GetAsMatrix(v_zg, n1, n2);

% Get the vector x
% A_{t} x = c_{t}
xk = SolveAx_b(Ak_fg,ck);


x = ...
    [
    xk(1:idx_col-1);
    0;
    xk(idx_col:end);
    ];



% Build the matrix Y_{t}
% Where Y(x) * z = E(z) * x
Yk = BuildY_RelativeDegree_STLN(x,m1,m2,n1,n2,k1,k2);

% Get vector of coefficients of f(x,y)
v_fxy = GetAsVector(fxy);

% Get vector of coefficients of g(x,y)
v_gxy = GetAsVector(gxy);

% Test Yk
test1 = Yk * [v_fxy;v_gxy];
test2 = Ak_fg * xk;
norm(test1-test2)

% Build the matrix P_{t}
% Where P * [f;g] = c_{t}
Pk = BuildP_RelativeDegree_STLN(m1,m2,n1,n2,idx_col,k1,k2);

% Test Pk
test1 = Pk * [v_fxy;v_gxy];
test2 = ck;
norm(test1-test2)


% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
xk = SolveAx_b(Ak_fg+Ak_zfzg,ck+hk);

% % Build the Matrix C consisting of H_{z} and H_{x}
H_z = Yk - Pk;

H_x = Ak_fg + Ak_zfzg;

C = [H_z H_x];

E = blkdiag( eye(nCoeffs_fxy + nCoeffs_gxy) , eye(nCoeffs_uxy + nCoeffs_vxy - 1));


% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
res_vec = (ck + hk) - (Ak_fg*xk);

start_point     =   ...
    [...
        z;...
        xk;
    ];

yy = start_point;

f = -(yy-start_point);

% Initialise the iteration counter
ite = 1;

% Set the termination criterion
condition(ite) = norm(res_vec)./ norm(ck);

while condition(ite) >  SETTINGS.MAX_ERROR_SNTLN &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN
    
    % Increment interation counter
    ite = ite + 1;
    
    % Get small petrubations by LSE problem ||Ey-f||  subject to  Cy=g
    y_lse = LSE(E,f,C,res_vec);
    
    % Increment cummulative peturbations
    yy = yy + y_lse;
    
    % obtain the small changes to z and x
    nEntries_z      = nCoeffs_fxy + nCoeffs_gxy;
    delta_zk        = y_lse(1:nEntries_z);
    delta_xk        = y_lse((nEntries_z+1):end);
    
    % Update z and x
    z = z + delta_zk;
    xk = xk + delta_xk;
    
    % Split vector z into vectors z_f and z_g
    v_zf = z(1 : nCoeffs_fxy);
    v_zg = z(nCoeffs_fxy + 1 : end);
    
    % Get zf as a matrix
    zf_xy = GetAsMatrix(v_zf, m1, m2);
    
    % Get zg as a matrix
    zg_xy = GetAsMatrix(v_zg, n1, n2);
    
    Bt = BuildT_Relative_Bivariate_2Polys(zf_xy, zg_xy, k1, k2);
    
    
    % Get the matrix E_{t} with optimal column removed
    Ak_zfzg = Bt;
    Ak_zfzg(:,idx_col) = [];
    
    % Get the column vector h_{t}, the optimal column removed from B_{t},
    % and equivalent to c_{t} removed from S_{t}
    hk = Bt(:,idx_col);
    
    % Get the updated vector x = [x1 x2] where x1 and x2 are vectors.
    % S(x1,x2)*[f;g] = ct
    % x_ls = SolveAx_b(At+Et,ct+ht);
    
    x = [...
        xk(1:idx_col-1);...
        0;...
        xk(idx_col:end)];
    
    % Build the matrix Y_{t} where Y_{t}(x)*z = E_{t}(z) * x
    Yk = BuildY_RelativeDegree_STLN(x, m1, m2, n1, n2, k1, k2);
    
    % Get the residual vector
    res_vec = (ck+hk) - ((Ak_fg + Ak_zfzg)*xk);
    
    % Update the matrix C
    H_z = Yk - Pk;
    H_x = Ak_fg + Ak_zfzg;
    C = [H_z H_x];
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    % Update the termination criterion
    condition(ite) = norm(res_vec)./norm(ck+hk) ;
    
    
end
SETTINGS.LOW_RANK_APPROX_REQ_ITE = ite;

LineBreakLarge();
fprintf([mfilename ' : ' sprintf('Required number of iterations: %i\n',ite)]);
LineBreakLarge();

PlotGraphs_STLN();

fxy_lr = fxy + zf_xy;
gxy_lr = gxy + zg_xy;

x = [...
        xk(1:idx_col-1);...
        -1;...
        xk(idx_col:end)];

vec_vxy = x(1:nCoeffs_vxy);
vec_uxy = -1.*x(nCoeffs_vxy+1:nCoeffs_vxy+nCoeffs_uxy);

uxy_lr = GetAsMatrix(vec_uxy,m1-k1,m2-k2);
vxy_lr = GetAsMatrix(vec_vxy,n1-k1,n2-k2);


end



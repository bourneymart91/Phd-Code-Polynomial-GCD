function [fxy_lr, gxy_lr, uxy_lr, vxy_lr] = STLN_Both(fxy, gxy, m, n, k, k1, k2, idx_col)
% STLN(fxy,gxy,m,n,t1,t2,opt_col)
%
% Given coefficients f(x,y) and g(x,y) find the low rank approximation of
% the Syvlester subresultant S_{t_{1},t_{2}}.
%
%
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y) 
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% k : (Int) Total Degree of d(x,y)
%
% k1 : (Int) Degree of d(x,y) with respect to x
%
% k2 : (Int) Degree of d(x,y) with respect to y
%
% idx_col : (Int) index of optimal column for removal from S_{t_{1},t_{2}}(f,g)
%
% % Outputs
%
% fxy_lr : (Matrix) Coefficients of f(x,y) with added perturbations
%
% gxy_lr : (Matrix) Coefficients of g(x,y) with added perturbations
%
% uxy_lr : (Matrix) Coefficients of u(x,y)
%
% vxy_lr : (Matrix) Coefficients of v(x,y)

global SETTINGS

% Get degree of polynomials f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degree of polynomial g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Get the number of coefficients in the polynomial f(x,y)
nCoefficients_fxy = (m1+1) * (m2+1);
nNonZeros_fxy = GetNumNonZeros(m1,m2,m);
nZeros_fxy = nCoefficients_fxy - nNonZeros_fxy;

% Get the number of coefficients in the polynomial g(x,y)
nCoefficients_gxy = (n1+1) * (n2+1);
nNonZeros_gxy = GetNumNonZeros(n1,n2,n);
nZeros_gxy = nCoefficients_gxy - nNonZeros_gxy;

% Get number of zeros in v
nCoefficients_vxy = (n1-k1+1) * (n2-k2+1);
nNonZeros_vxy = GetNumNonZeros(n1-k1,n2-k2,n-k);
nZeros_vxy = nCoefficients_vxy - nNonZeros_vxy;


% Get number of zeros in u(x,y)
nCoefficients_uxy = (m1-k1+1) * (m2-k2+1);
nNonZeros_uxy = GetNumNonZeros(m1-k1, m2-k2, m-k);
nZeros_uxy = nCoefficients_uxy - nNonZeros_uxy;

% Build the matrix T_{n1-k1,n2-k2}(f)
Tf = BuildT1_Both_Bivariate(fxy, m, n-k, n1-k1, n2-k2);

% Build the matrix T_{m1-k1,m2-k2}(g)
Tg = BuildT1_Both_Bivariate(gxy, n, m-k, m1-k1, m2-k2);

% Remove the columns of T(f) and T(g) which correspond to the zeros in u(x,y)
% and v(x,y) which are removed from the solution vector x.

% %
% %
% Remove the extra rows of T1 and T2 associated with zeros of f*v and g*u
% Build the matrix S_{t,t_{1},t_{2}} with reduced rows and columns
Sk_fg = [Tf Tg];

% Remove optimal column
Ak_fg = Sk_fg;
Ak_fg(:,idx_col) = [];
ck = Sk_fg(:,idx_col);

% Build the matrix Et
Ak_zfzg = zeros(size(Ak_fg));
hk = zeros(size(ck));

% Build the vector z, which is the vector of perturbations of f and g.
z = zeros(nNonZeros_fxy + nNonZeros_gxy,1);

% Get the vector of coefficients zf
v_zf = z(1:nNonZeros_fxy);

% Get the vector of coefficeints zg
v_zg = z(nNonZeros_fxy + 1:end);

% Get zf as a matrix
% EDIT 10/03/2016 - vZ_fxy has zeros removed, so include the zeros to form
% a matrix m1+1 * m2+1
mat_zf = GetAsMatrix_Version1(...
    [
    v_zf;
    zeros(nZeros_fxy,1)
    ],m1,m2);

% Get zg as a matrix
mat_zg = GetAsMatrix_Version1(...
    [
    v_zg;
    zeros(nZeros_gxy,1)
    ],n1,n2);

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
Yk = BuildY_BothDegree_STLN(x, m, m1, m2, n, n1, n2, k, k1, k2);

% Test Y_{k}
v_fxy = GetAsVector_Version1(fxy);
v_fxy = v_fxy(1:nNonZeros_fxy,:);

v_gxy = GetAsVector_Version1(gxy);
v_gxy = v_gxy(1:nNonZeros_gxy,:);

test1a = Yk * [v_fxy;v_gxy];
test1b = Sk_fg * x;
test1 = norm(test1a-test1b);
display(test1)

% Build the matrix P_{t}
% Where P * [f;g] = c_{t}
Pk = BuildP_BothDegree_STLN(m, m1, m2, n, n1, n2, k, k1, k2, idx_col);
test2a = Pk * [v_fxy;v_gxy];
test2b = ck;
test2 = norm(test2a-test2b);
display(test2)

% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
xk = SolveAx_b(Ak_fg+Ak_zfzg,ck+hk);

H_z = Yk - Pk;

H_x = Ak_fg + Ak_zfzg;

C = [H_z H_x];

E = blkdiag( eye(nNonZeros_fxy + nNonZeros_gxy) , eye(nNonZeros_uxy + nNonZeros_vxy - 1));


% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
res_vec = (ck + hk) - (Ak_fg*xk);

start_point     =   ...
    [...
    z;...
    xk;
    ];

yy = start_point;

f = -(yy - start_point);

% Initialise the iteration counter
ite = 1;

% Set the termination criterion
condition(ite) = norm(res_vec)./ norm(ck);

while condition(ite) >  SETTINGS.MAX_ERROR_SNTLN &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN
    
    % Increment interation counter
    ite = ite + 1;
    
    % Get small petrubations by LSE problem ||Ey-f||  subject to  Cy=g
    y_lse = LSE(E, f, C, res_vec);
    
    % Increment cummulative peturbations
    yy = yy + y_lse;
    
    % obtain the small changes to z and x
    nEntries_z      = nNonZeros_fxy + nNonZeros_gxy;
    delta_zk        = y_lse(1 : nEntries_z);
    delta_xk        = y_lse((nEntries_z+1):end);
    
    % Update z and x
    z     = z + delta_zk;
    xk    = xk + delta_xk;
    
    % Split vector z into vectors z_f and z_g
    v_zf = z(1 : nNonZeros_fxy);
    v_zg = z(nNonZeros_fxy + 1:end);
    
    % Get zf as a matrix
    mat_zf = GetAsMatrix_Version1(...
        [
        v_zf;
        zeros(nZeros_fxy,1)
        ]...
        ,m1,m2);
    
    % Get zg as a matrix
    mat_zg = GetAsMatrix_Version1(...
        [
        v_zg;
        zeros(nZeros_gxy,1)
        ]...
        ,n1,n2);
    
    % Build the matrix B_{t} = [E1(zf) E2(zg)]
    T1_zf = BuildT1_Both_Bivariate(mat_zf, m, n-k, n1-k1, n2-k2);
    T1_zg = BuildT1_Both_Bivariate(mat_zg, n, m-k, m1-k1, m2-k2);
    
    % Build the matrix B_{t} equivalent to S_{t}
    St_zfzg = [T1_zf T1_zg];
    
    % Get the matrix E_{t} with optimal column removed
    Ak_zfzg = St_zfzg;
    Ak_zfzg(:,idx_col) = [];
    
    % Get the column vector h_{t}, the optimal column removed from B_{t},
    % and equivalent to c_{t} removed from S_{t}
    hk = St_zfzg(:,idx_col);
    
    
    x = [...
        xk(1:idx_col-1);...
        0;...
        xk(idx_col:end)];
    
    % Build the matrix Y_{t} where Y_{t}(x)*z = E_{t}(z) * x
    Yk = BuildY_BothDegree_STLN(x, m, m1, m2, n, n1, n2, k, k1, k2);
    
    % Get the residual vector
    res_vec = (ck + hk) - ((Ak_fg + Ak_zfzg) * xk);
    
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
LineBreakLarge()
fprintf([mfilename ' : ' sprintf('\nRequired number of iterations: %i\n',ite)]);
LineBreakLarge()

fxy_lr = fxy + mat_zf;
gxy_lr = gxy + mat_zg;

x = [...
    xk(1:idx_col-1);...
    -1;...
    xk(idx_col:end)...
    ];

vec_vxy = x(1:nNonZeros_vxy);
vec_uxy = -1.* x(nNonZeros_vxy+1:end);

vxy_lr = GetAsMatrix_Version1([vec_vxy ; zeros(nZeros_vxy,1)], n1-k1, n2-k2);
uxy_lr = GetAsMatrix_Version1([vec_uxy ; zeros(nZeros_uxy,1)], m1-k1, m2-k2);

PlotGraphs_STLN()


end


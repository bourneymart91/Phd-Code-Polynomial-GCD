function [fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
    SNTLN(fxy, gxy, i_alpha, i_th1, i_th2, k)
% SNTLN(fxy,gxy,alpha,th1,th2,k)
%
% Compute the low rank approximation of the Sylvester matrix S_{k}(f,g)
% by method of Structured Total Least Norm STLN.
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) in the Bernstein basis
%
% gxy : (Matrix) Coefficients of polynomial g(x,y) in the Bernstein basis
%
% alpha : (Float) Optimal value of \alpha
%
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}
%
% k : (Int) Degree of common divisor d(x,y) and index of the Sylvester matrix
%     S_{k} whose low rank approximation is considered.
%
% % Outputs.
%
% fxy_lr : (Matrix) Coefficients f(x,y) obtained from the low rank approx of S_{k}(f,g)
%
% gxy_lr : (Matrix) Coefficients g(x,y) obtained from the low rank approx of S_{k}(f,g)
%
% uxy_lr : (Matrix) Coefficients u(x,y) obtained from the low rank approx of S_{k}(f,g)
%
% vxy_lr : (Matrix) Coefficients v(x,y) obtained from the low rank approx of S_{k}(f,g)
%
% alpha_lr : (Float)
%
% th1_lr : (Float)
%
% th2_lr : (Float)

global SETTINGS

% Set the initial iteration number
ite = 1;

% Set initial values of \alpha, \theta_{1} and \theta_{2}
th1(ite) = i_th1;
th2(ite) = i_th2;
alpha(ite) = i_alpha;

% Get the degree of f(x,y) and g(x,y)
[m,~] = GetDegree_Bivariate(fxy);
[n,~] = GetDegree_Bivariate(gxy);

% Get number of coefficients in f(x,y)
nCoefficients_fxy = nchoosek(m+2,2);

% Get the number of zeros in the matrix containing f(x,y)
nZeros_fxy = nchoosek(m+1,2);

% Get number of coefficients in g(x,y)
nCoefficients_gxy = nchoosek(n+2,2);

% Get the number of zeros in the matrix containing g(x,y)
nZeros_gxy = nchoosek(n+1,2);

% Get the number of coefficients in both f(x,y) and g(x,y)
nCoefficients_fg = nCoefficients_fxy + nCoefficients_gxy;

% Get the number of coefficients in v(x,y)
nCoefficients_vxy = nchoosek(n-k+2,2);

% Get the number of coefficients in u(x,y)
nCoefficients_uxy = nchoosek(m-k+2,2);

% Get the number of coefficients in the unknown vector x, where A_{k}x = c_{k}
nCoefficients_x = nCoefficients_uxy + nCoefficients_vxy - 1;

% Get the number of columns in C_{n-k}(f), the first partition of the
% Sylvester subresultant matrix S_{k}(f,g)
nColumns_Tf = nCoefficients_vxy;

% Get the number of columns in C_{m-k}(f), the second partititon of the
% Sylvester subresultant matrix S_{k}(f,g)
nColumns_Tg = nCoefficients_uxy;

% Get the total number of columns in the Sylvester subresultant matrix
% S_{k}(f,g)
nColumns_Sk = nColumns_Tf + nColumns_Tg;


% %
% Preprocessing

% Get the polynomial f(\omega_{1},\omega_{2})
fww = GetWithThetas(fxy, m, th1(ite), th2(ite));

% Get the polynomial g(\omega_{1},\omega_{2})
gww = GetWithThetas(gxy, n, th1(ite), th2(ite));

% % 
% Build matrices

% Build the kth Sylvester Subresultant matrix S_{k}(f,g) = D*T(f,g)*Q

% Build the diagonal matrix D^{-1}
D = BuildD_2Polys(m, n - k);

% Build the matrix T_{n-k}(f)
T1_fww = BuildT1(fww, m, n - k);

% Build the matrix T_{m-k}(g)
T1_gww = BuildT1(gww, n, m - k);

% Build the matrix Q_{k}
Q = BuildQ_2Polys(m, n, k);

% Build the kth Sylvester subresultant matrix.
DTQ_fg = D*[T1_fww alpha.*T1_gww] * Q;

% Get index of optimal column for removal from S_{k}(f,g)
idx_col = GetOptimalColumn(DTQ_fg);

% Get c_{k} optimal column removed from S_{k}(f,g)
ck = DTQ_fg(:,idx_col);

% Create the identity matrix
I = eye(nColumns_Sk, nColumns_Sk);
M = I;
M(:,idx_col) = [];
e = I(:,idx_col);


% %  

% Get partial derivatives of f(\omega_{1},\omega_{2}) with respect to \alph
fww_wrt_alpha = zeros(m + 1, m + 1);

% Get partial derivatives of g(\omega_{1},\omega_{2}) with respect to \alph
alpha_gww_wrt_alpha = gxy;

% Calculate the partial derivatives of f(w,w) with respect to theta_1
fww_wrt_th1 = Differentiate_wrt_th1(fww, th1(ite));

% Calculate the partial derivative of g(w,w) with respect to theta_1
gww_wrt_th1 = Differentiate_wrt_th1(gww, th1(ite));

% Calculate the partial derivative of f(w,w) with respect to theta_2
fww_wrt_th2 = Differentiate_wrt_th2(fww, th2(ite));

% Calculate the partial deriviates of g(w,w) with respect to theta_2
gww_wrt_th2 = Differentiate_wrt_th2(gww, th2(ite));

% Build the derivative of T(f,g) with respect to alpha
DTQ_alpha = BuildDTQ_2Polys(fww_wrt_alpha, alpha_gww_wrt_alpha, k);

% Calculate the derivative of T(f,g) with respect to theta_{1}
DTQ_wrt_th1 = BuildDTQ_2Polys(fww_wrt_th1, alpha(ite).*gww_wrt_th1, k);

% Calcualte the derivative of T(f,g) with respect to theta_2
DTQ_wrt_th2 = BuildDTQ_2Polys(fww_wrt_th2, alpha(ite).*gww_wrt_th2, k);


% 
% Initialise the vector z of structured perturbations
% if we are working with strictly the roots problem, the number of entries
% in z can be reduced.

nRows_Sk_fg = nchoosek(m + n - k + 2, 2);

zk = zeros(nCoefficients_fg , 1);

% 
% Initilaise the derivative of N wrt alpha.
DNQ_wrt_alpha   = zeros(nRows_Sk_fg, nColumns_Tf + nColumns_Tg);

% Initilaise the derivative of N wrt theta_1.
DNQ_wrt_th1   = zeros(nRows_Sk_fg, nColumns_Tf + nColumns_Tg);

% Initialise the derivative of N wrt theta 2
DNQ_wrt_th2   = zeros(nRows_Sk_fg, nColumns_Tf + nColumns_Tg);

% %
% Initialise the derivative of h_{k}
% Calculate the derivatives wrt alpha and theta of the column of DNQ
% that is moved to the right hand side.
hk_wrt_alpha = DNQ_wrt_alpha*e;
hk_wrt_th1 = DNQ_wrt_th1*e;
hk_wrt_th2 = DNQ_wrt_th2*e;

% %
% Build the matrix DPG
DPG = BuildDPG_SNTLN(m,n,k,alpha(ite),th1(ite),th2(ite),idx_col);
f = GetAsVector(fxy);
f = f(1:nCoefficients_fxy);

g = GetAsVector(gxy);
g = g(1:nCoefficients_gxy);

test1a = DPG*[f;g];
test1b = ck;
test1 = (test1a - test1b);
display(norm(test1));

% %
% Calculate the derivatives wrt alpha and theta of the removed column.

% Get parital c_{k} with respect to \alpha 
ck_wrt_alpha = DTQ_alpha(:,idx_col);

% Get partial c_{k} with respect to \theta_{1}
ck_wrt_th1 = DTQ_wrt_th1(:,idx_col);

% Get partial c_{k} with respect to \theta_{2}
ck_wrt_th2 = DTQ_wrt_th2(:,idx_col);

% Remove optimal column from S_{k}(f,g)
Ak_fg = DTQ_fg;
Ak_fg(:,idx_col) = [];



% Get least squares solution x_ls
xk = SolveAx_b(Ak_fg,ck);

% Obtain the vector \hat{x} which contains x_ls with a zero inserted into
% the position of the optimal column.
x = [xk(1:idx_col-1) ; 0 ; xk(idx_col:end)];

% Build Matrix Y, where Y(v,u)*[f;g] = S(f,g)*[u;v]
DYQ = BuildDYQ_SNTLN(m, n, k, x, alpha(ite), th1(ite), th2(ite));
test2a = DYQ * [f;g];
test2b = DTQ_fg * [x];
test2 = (test2a - test2b);
display(norm(test2));


% Calculate the initial residual r = ck - (Ak*x)
res_vec = ck - (DTQ_fg*M*xk);

% % Get the matrix p, which will store all the perturbations returned from LSE file
nEntries = nCoefficients_fxy + nCoefficients_gxy + nCoefficients_x + 3;

% Set the intial value of E to the identity matrix
E = eye(nEntries);

% Create the matrix D(T+N)Q, initially N is empty so this is the same as T.
DTNQ = BuildDTQ_2Polys(fww, alpha(ite).*gww,k);

% Create The matrix (T+N) with respect to alpha
DTNQ_wrt_alpha = BuildDTQ_2Polys(fww_wrt_alpha, alpha_gww_wrt_alpha,k);

% Create The matrix (T+N) with respect to \theta_{1}
DTNQ_wrt_th1 = BuildDTQ_2Polys(fww_wrt_th1,alpha(ite).*gww_wrt_th1,k);

% Create The matrix (T+N) with respect to \theta_{2}
DTNQ_wrt_th2 = BuildDTQ_2Polys(fww_wrt_th2,alpha(ite).*gww_wrt_th2,k);

% %
% Create the matrix C for input into iteration

H_z = DYQ - DPG;

H_x = DTNQ*M;

H_alpha  = DTNQ_wrt_alpha*M*xk - ...
    (ck_wrt_alpha + hk_wrt_alpha);

H_th1 = DTNQ_wrt_th1*M*xk - ...
    (ck_wrt_th1 + hk_wrt_th1);

H_th2 = DTNQ_wrt_th2*M*xk - ...
    (ck_wrt_th2 + hk_wrt_th2);

C = [H_z H_x H_alpha H_th1 H_th2];

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    zk;...
    xk;...
    alpha(ite);...
    th1(ite);...
    th2(ite)
    ];

yy  = start_point;

f = -(yy - start_point);

% Set the termination criterion to a large value. It will be
% over written later.
vCondition(ite) = norm(res_vec)/norm(ck);


while vCondition(ite) > SETTINGS.STLN_MAX_ERROR && ite < SETTINGS.STLN_MAX_ITERATIONS
    
    
    % Get small petrubations by LSE
    y = LSE(E,f,C,res_vec);
    
    % Increment interation counter
    ite = ite + 1;
    
    % Increment cummulative peturbations
    yy = yy + y;
    
    % Get the entries corresponding to perturbations of f(x,y) and g(x,y)
    delta_zk = y(1 : nCoefficients_fxy + nCoefficients_gxy);
    
    % Remove the zk coefficients from the vector y.
    y(1:nCoefficients_fxy + nCoefficients_gxy) = [];
    
    % Get the entries of y corresponding to x
    delta_xk = y(1:nCoefficients_x,1);
    
    % Remove the coefficients from y
    y(1:nCoefficients_x) = [];
    
    % Get the entry corresponding to \alpha
    delta_alpha = y(1:1);
    
    % Remove the entry from the list
    y(1) = [];
    
    % Get the entry in y corresponding to \theta_{1}
    delta_th1 = y(1:1);
    
    % Remove the entry from vector y
    y(1) = [];
    
    % Get the entry in y corresponding to \theta_{2}
    delta_th2 = y(1:1);
    
    % Remove the entry from vector y
    y(1) = [];
    
    % % Update the variables
    
    % Update variables z_{k}, where z_{k}, where z_{k} are perturbations in
    % the coefficients of f and g.
    zk = zk + delta_zk;
    
    % Update x_{k}
    xk = xk + delta_xk;
    
    % Update alpha
    alpha(ite) = alpha(ite-1) + delta_alpha;
    
    % Update \theta_{1}
    th1(ite) = th1(ite-1) + delta_th1;
    
    % Update \theta_{2}
    th2(ite) = th2(ite-1) + delta_th2;
    
    % %
    % Obtain polynomials in modified Bernstein basis a_{i}\theta^{i}
    
    % Get new f(\omega_{1},\omega_{2}) with improved \theta_{1} and
    % \theta_{2}.
    fww = GetWithThetas(fxy, m, th1(ite), th2(ite));
    
    % Get new g(w,w)
    gww = GetWithThetas(gxy, n, th1(ite), th2(ite));
    
    % Construct the kth Sylvester subresultant matrix DTQ.
    DTQ_fg = BuildDTQ_2Polys(fww, alpha(ite).*gww,k);
    
    ck = DTQ_fg*e;
    
    % %
    % Get partial derivatives
    
    % Get the partial derivative of f(\omega_{1},\omega_{2}) with respect
    % to \alpha
    fww_wrt_alpha = zeros(m+1,m+1);
    
    % Get the partial derivative of g(\omega_{1},\omega_{2}) with respect
    % to \alpha
    alpha_gww_wrt_alpha = gww;
    
    % Calculate the partial derivatives of f(\omega_{1},\omega_{2}) with
    % respect to \theta_{1}
    
    % Get the partial derivative of f(w1,w2) with respect to theta_{1}
    fww_wrt_th1 = Differentiate_wrt_th1(fww,th1(ite));
    
    % Get the partial derivative of g(w1,w2) with respect to theta_{1}
    gww_wrt_th1 = Differentiate_wrt_th1(gww,th1(ite));
    
    % Get the partial derivative of f with respect to theta_{2}
    fww_wrt_th2 = Differentiate_wrt_th2(fww,th2(ite));
    
    % Get the partial derivative of g with respect to theta_{2}
    gww_wrt_th2 = Differentiate_wrt_th2(gww,th2(ite));
    
    % Calculate the Partial derivative of T with respect to alpha.
    DTQ_wrt_alpha = BuildDTQ_2Polys(fww_wrt_alpha, alpha_gww_wrt_alpha,k);
    
    % Calculate the partial derivative of DTQ with respect to theta_{1}
    DTQ_wrt_th1 = BuildDTQ_2Polys(fww_wrt_th1,alpha(ite).*gww_wrt_th1,k);
    
    % Calculate the partial derivative of DTQ with respect to theta_{2}
    DTQ_wrt_th2 = BuildDTQ_2Polys(fww_wrt_th2,alpha(ite).*gww_wrt_th2,k);

    % Calculate the derivatives of c_{k} with respect to \alpha, \theta_{1}
    % and \theta_{2}
    ck_wrt_alpha = DTQ_wrt_alpha*e;
    ck_wrt_th1 = DTQ_wrt_th1*e;
    ck_wrt_th2 = DTQ_wrt_th2*e;
    
    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    vec_z_fxy      = zk(1:nCoefficients_fxy);
    vec_z_gxy      = zk(nCoefficients_fxy + 1 :end);
    
    % Get the vectors z_fx and z_gx as matrices, which match the shape of
    % f(x) and g(x).
    z_fxy = GetAsMatrix([vec_z_fxy ; zeros(nZeros_fxy,1)], m, m);
    z_gxy = GetAsMatrix([vec_z_gxy ; zeros(nZeros_gxy,1)], n, n);
    
    % Get matrices z_fw_mat and z_gw_mat, by multiplying rows by
    % theta_{1}^{i} and columns by theta_{2}^{j}
    z_fww = GetWithThetas(z_fxy, m, th1(ite), th2(ite));
    z_gww = GetWithThetas(z_gxy, n, th1(ite), th2(ite));
    
    % Calculate the derivatives of z_fw and z_gw with repect to \alpha.
    zfw_wrt_alpha = zeros(m+1,m+1);
    alpha_zgw_wrt_alpha = z_gww;
    
    % Calculate the derivative of z_fw with respect to \theta_{1}.
    zfw_wrt_th1 = Differentiate_wrt_th1(z_fww, th1(ite));
    
    % Calculate the derivative of z_fw with respect to \theta_{2}
    zfw_wrt_th2 = Differentiate_wrt_th2(z_fww, th2(ite));
    
    % Calculate the derivative of z_gw with respect ot theta1
    zgw_wrt_th1 = Differentiate_wrt_th1(z_gww, th1(ite));
    
    % Calculate the deriviate of z_gw with respect to theta2
    zgw_wrt_th2 = Differentiate_wrt_th2(z_gww, th2(ite));
    
    % Build the coefficient Matrix N = [T(z1) T(z2)], of structured perturbations, with
    % same structure as DTQ.
    DNQ = BuildDTQ_2Polys(z_fww,alpha(ite).*z_gww,k);
    
    % Build the coefficient matrix N with respect to alpha
    DNQ_wrt_alpha = BuildDTQ_2Polys(zfw_wrt_alpha, alpha_zgw_wrt_alpha,k);
    
    % Calculate the derivatives of DNQ with respect to theta
    DNQ_wrt_th1 = BuildDTQ_2Polys(zfw_wrt_th1, alpha(ite).*zgw_wrt_th1,k);
    
    % Calculate the derivatives of DNQ with respect to theta
    DNQ_wrt_th2 = BuildDTQ_2Polys(zfw_wrt_th2, alpha(ite).*zgw_wrt_th2,k);
    
    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{k} the column of S_{k} moved to the RHS
    hk = DNQ*e;
    
    % Calculate the derivative of h with respect to \alpha
    hk_alpha = DNQ_wrt_alpha*e;
    
    % Calculate the derivative of h with respect to \theta_{1}
    hk_th1 = DNQ_wrt_th1*e;
    
    % Calculate the derivative of h with respect to \theta_{2}
    hk_th2 = DNQ_wrt_th2*e;
    
    % Build the matrix D*(T+N)*Q
    DTNQ = BuildDTQ_2Polys(fww + z_fww, alpha(ite).*(gww + z_gww), k);
    
    % Calculate the paritial derivative of D*(T+N)*Q with respect to
    % \alpha
    DTNQ_alpha = BuildDTQ_2Polys(fww_wrt_alpha + zfw_wrt_alpha,...
        alpha_gww_wrt_alpha + alpha_zgw_wrt_alpha,...
        k);
    
    
    % Calculate the paritial derivative of D*(T+N)*Q with respect to
    % \theta_{1}
    DTNQ_th1 = BuildDTQ_2Polys(...
        fww_wrt_th1 + zfw_wrt_th1,...
        alpha(ite).*(gww_wrt_th1 + zgw_wrt_th1),...
        k);
    
    % Calculate the paritial derivative of D*(T+N)*Q with respect to
    % \theta_{2}
    DTNQ_th2 = BuildDTQ_2Polys(...
        fww_wrt_th2 + zfw_wrt_th2,...
        alpha(ite).*(gww_wrt_th2 + zgw_wrt_th2),...
        k);
    
    % Calculate the matrix DYQ where DYQ is the matrix such that 
    % D*Y(x1,x2)*Q*[f;g] = D*S(f,g)*Q*[x1;x2]

    % Insert a zero into the position of the optimal_column.
    % Partition x_ls into the two parts for coefficients of
    first_part = xk(1:(idx_col-1));
    second_part = xk(idx_col:end);
    x = [first_part ; 0 ; second_part];
    
    % Build the matrix DYQ
    DYQ = BuildDYQ_SNTLN(m, n, k, x, alpha(ite), th1(ite), th2(ite));

    % Calculate the matrix DPG
    DPG = BuildDPG_SNTLN(m, n, k, alpha(ite), th1(ite), th2(ite), idx_col);
       
    % Get residual as a vector
    res_vec = (ck+hk) - DTNQ*M*xk ;
    
    % %
    % Create the matrix C. This is made up of five submatrices, H_{z}, H_{x},
    % H_{\alpha}, H_{\theta_{1}} and H_{\theta_{2}}
    
    Hz = DYQ - DPG;
    
    Hx = DTNQ*M;
    
    H_alpha = DTNQ_alpha*M*xk - (ck_wrt_alpha + hk_alpha);
    
    H_th1 = DTNQ_th1*M*xk - (ck_wrt_th1 + hk_th1);
    
    H_th2 = DTNQ_th2*M*xk - (ck_wrt_th2 + hk_th2);
    
    C = [Hz,Hx,H_alpha,H_th1, H_th2];  % the matrix C
           
    % Calculate the normalised residual of the solution.
    vCondition(ite) = norm(res_vec) / norm(ck + hk);
    
    % Update fnew - used in LSE Problem.
    f = -(yy - start_point);
    
end

%
plotgraphs_STLN();

%
LineBreakLarge()
fprintf([mfilename ' : ' sprintf('SNTLN Number of iterations required : %i \n',ite)])
LineBreakLarge()
SETTINGS.LOW_RANK_APPROX_REQ_ITE = ite;

% Update z and x_ls
vec_z = yy(1 : nCoefficients_fxy + nCoefficients_gxy);

% Get z1 and z2
zf = vec_z(1:nCoefficients_fxy);
zg = vec_z(nCoefficients_fxy+1:end);

% Get Coefficients of polynomials f(x,y) and g(x,y)
vec_zf = [zf ; zeros(nZeros_fxy,1)];
vec_zg = [zg ; zeros(nZeros_gxy,1)];

zf_xy = GetAsMatrix(vec_zf, m, m);
zg_xy = GetAsMatrix(vec_zg, n, n);


fxy_lr = fxy + zf_xy;
gxy_lr = gxy + zg_xy;

% Get coefficients of polynomials xu(x,y) and xv(x,y)
% split the vector x into \hat{x}_{1} and \hat{x}_{2}
% Get number of coefficients in x1

% split the vector x into \hat{x}_{1} and \hat{x}_{2}
% Get number of coefficients in x1
nCoeffs_xv = nchoosek(n-k+2, 2);
nCoeffs_xu = nchoosek(m-k+2, 2);

vec_xvxu = [xk(1:idx_col-1) ; -1 ; xk(idx_col:end)];

xv = vec_xvxu(1:nCoeffs_xv);
xu = -1.*vec_xvxu(nCoeffs_xv + 1 : nCoeffs_xv + nCoeffs_xu);

% Get x1 as a matrix of coefficients for input into BuildT1() function

nZeros_xv = nchoosek(n-k+1, 2);
nZeros_xu = nchoosek(m-k+1, 2);

% Get vectors of coefficients of v and u
vec_xv = [ xv ; zeros(nZeros_xv, 1)];
vec_xu = [ xu ; zeros(nZeros_xu, 1)];

% Get as matrix
mat_xv = GetAsMatrix(vec_xv, n-k, n-k);
mat_xu = GetAsMatrix(vec_xu, m-k, m-k);

vxy_lr = mat_xv;
uxy_lr = mat_xu;

% %
%
alpha_lr = alpha(ite);
th1_lr = th1(ite);
th2_lr = th2(ite);





end




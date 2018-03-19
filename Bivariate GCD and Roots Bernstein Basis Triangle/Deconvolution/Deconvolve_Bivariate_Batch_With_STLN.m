function [arr_hxy] = Deconvolve_Bivariate_Batch_With_STLN(arr_fxy, vDeg_fxy)
% Get the set of polynomials h_{i} given by the deconvolution of the
% polynomials f_{i}, where h_{i} = f_{i-1}/f_{i}
%
% % Inputs.
%
%
% arr_fxy : Array of polynomials f_{i}(x,y)
%
% vDeg_fxy : Vector containing the total degree of the polynomials f{i}
%
%
% % Outputs.
%
%
% arr_hxy :

% Global Variables
global SETTINGS


% %
% %
% Form the left hand side matrix

% Get number of polynomials in the array arr_fxy of f_{i}(x,y)
% Get number of polynomials in the array arr_hxy of h_{i}(x,y)
nPolys_arr_fxy = size(arr_fxy, 1);
nPolys_arr_hxy = nPolys_arr_fxy - 1;

% %
% %
% Get the degree of each polynomial f_{i}(x,y)

% Initialise vector to store degrees of f_{i}(x,y)
vDeg_arr_fxy = zeros(nPolys_arr_fxy,1);

% For each polynomial f_{i} get its degree
for i = 1:1:nPolys_arr_fxy
    
    vDeg_arr_fxy(i)  = GetDegree_Bivariate(arr_fxy{i});
    
end

% %
% %
% Get the degrees of each polynomial h_{i}(x,y) = f_{i-1}/f_{i}
vDeg_arr_hxy = vDeg_arr_fxy(1:end-1) - vDeg_arr_fxy(2:end);

% %
% %
% Preprocess polynomials f_{i}(x,y) based on their occurences in the LHS
% matrix
%
% y - Preprocess.
% n - Dont preprocess.

if( SETTINGS.PREPROC_DECONVOLUTIONS)
    
    [th1, th2] = GetOptimalTheta(arr_fxy);
    
else
    th1 = 1;
    th2 = 1;
    
end

% Get preprocessed f(\omega_{1},\omega_{2}) from f(x,y)

arr_fww = GetPolynomialArrayWithThetas(arr_fxy, th1, th2);


% %
% %
% %
% Write deconvolutions in form D^{-1}C(f)Q] h = RHS_f


% Form the right hand side vector
vRHS_fww = BuildRHS_vec(arr_fww);

% Build the matrix C(f1,...,fd)

DT_fwwQ = BuildLHS_Matrix(arr_fww);

% Get vector of coefficients of the polynomials h_{i}(x,y)
v_hww = SolveAx_b(DT_fwwQ,vRHS_fww);

% Get array of polynomials
arr_hww = GetPolynomialArrayFromVector(v_hww,vDeg_arr_hxy);


% Let z be vector of perturbations of polynomails f_{i} such that
% z = [z0 z1 ... zd]
arr_zww = cell(nPolys_arr_fxy,1);

for i = 1 : 1 : nPolys_arr_fxy
    m = vDeg_arr_fxy(i);
    arr_zww{i,1} = zeros(m+1,m+1);
end

% Build vector z, consisting of all vectors z_{i}
v_zww = Get_vec_z(arr_zww);

v_RHS_zww = BuildRHS_vec(arr_zww);

% %
% %
% Build the matrix P

% Get number of coefficients across all polynomials f_{i}(x,y)
v_nCoefficients_arr_fxy = zeros(nPolys_arr_fxy,1);

for i = 1:1:nPolys_arr_fxy
    m = vDeg_arr_fxy(i);
    v_nCoefficients_arr_fxy(i,1) = nchoosek(m+2,2);
end

% Get number of coefficients across all polynomials h_{i}(x,y)
v_nCoefficients_arr_hxy = zeros(nPolys_arr_hxy,1);

for i = 1:1:nPolys_arr_hxy
    n = vDeg_arr_hxy(i);
    v_nCoefficients_arr_hxy(i,1) = nchoosek(n+2,2);
end

% Get number of coefficients in all polynomials f_{i}(x,y)
nCoefficients_fxy = sum(v_nCoefficients_arr_fxy);
nCoefficients_hxy = sum(v_nCoefficients_arr_hxy);

% Get sum of all coefficients except the final f_{i}(x,y)
nCoefficients_rhs = sum(v_nCoefficients_arr_fxy(1:end-1));

% Finally build the matrix P
P = [...
    eye(nCoefficients_rhs) ...
    zeros(nCoefficients_rhs,v_nCoefficients_arr_fxy(end))
    ];

% %
% %
% %
% Build matrix Y where E(z)*h = Y(h)*z
DY_hQ = BuildDYQ(arr_hww,vDeg_arr_fxy);

% Compute the initial residual
t = ((vRHS_fww + v_RHS_zww) - (DT_fwwQ*v_hww));

% Set the iteration number
ite = 1;

% Define the weight vector
E = [eye(nCoefficients_fxy) zeros(nCoefficients_fxy, nCoefficients_hxy)];

% Build the matrix C for the LSE problem
C_z = (DY_hQ) - P;
C_h = DT_fwwQ;
C = [ C_z C_h ];


% Get the condition 
condition(ite) = norm(t) ./ norm(vRHS_fww);

start_point = ...
    [
        v_zww;
        v_hww;
    ];

yy = start_point;



s = E * (start_point - yy);

% Perform iteration to obtain perturbations

while (condition(ite) > SETTINGS.MAX_ERROR_DECONVOLUTIONS)  && ...
        (ite < SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS)
    
    % Use the QR decomposition to solve the LSE problem and then
    % update the solution.
    % min |Fy-s| subject to Gy=t
    y = LSE_new(E, s, C, t);
    
    yy = yy + y;
    
    % Output y gives delta h and delta z
    delta_z = y(1 : nCoefficients_fxy);
    delta_h = y(nCoefficients_fxy + 1 : end);
    

    % Add structured perturbations to vector z.
    v_zww = v_zww + delta_z;
    
    % Add structured perturbations to vector h.
    v_hww = v_hww + delta_h;
    
    arr_zww = GetPolynomialArrayFromVector(v_zww, vDeg_arr_fxy);
    arr_hww = GetPolynomialArrayFromVector(v_hww, vDeg_arr_hxy);
    
    
    % Increment s in LSE Problem
    s = E * (start_point - yy);
    
    % Build the iterative DY_hQ
    DY_hQ = BuildDYQ(arr_hww, vDeg_arr_fxy);
    
    % Build D[C(f)+C(h)]Q
    DC_fQ = BuildLHS_Matrix(arr_fww);
    DC_zQ = BuildLHS_Matrix(arr_zww);
    
    % Build G
    C = [(DC_fQ + DC_zQ) (DY_hQ - P)];
    
    % Update the RHS_vector
    vRHS_fww = BuildRHS_vec(arr_fww);
    vRHS_zww = BuildRHS_vec(arr_zww);
    
    % Calculate residual and increment t in LSE Problem
    t = ((vRHS_fww + vRHS_zww ) - ((DC_fQ + DC_zQ)*v_hww));
    
    
    % Get the condition
    condition(ite +1) = norm(t)./norm((vRHS_fww + vRHS_zww));
    
    % Increment iteration number
    ite = ite + 1;
    
end

% Get array of polynomials h_{i}(x,y) from h_{i}(\omega1,\omega2)
arr_hxy = GetPolynomialArrayWithoutThetas(arr_hww,th1,th2);

LineBreakLarge()
fprintf([mfilename ': Iterations Required : ' num2str(ite) '\n'])
LineBreakLarge()


% Plot Residuals
figure_name = sprintf([mfilename ' : Residuals']);
figure('Name',figure_name)
hold on
plot(log10(condition),'-s','LineWidth',2)
xlabel('$ite$','Interpreter','latex','FontSize', 20)
ylabel('$\log_{10}$','Interpreter','latex','FontSize', 20)
hold off


end

function C_fxy = BuildLHS_Matrix(arr_fxy)
% Build the matrix C(f1,...,fd)
%
% Inputs.
%
% arr_fxy : Array of polynomials f(x,y)
%
% Outputs.
%
% C_fxy :

% Get number of polynomials in f_{i}(x,y)
nPolys_arr_fxy = size(arr_fxy, 1);

% Get degree of each polynomial f_{i}(x,y)
vDeg_arr_fxy = zeros(nPolys_arr_fxy, 1);

for i = 1 : 1 : nPolys_arr_fxy
    
    vDeg_arr_fxy(i) = GetDegree_Bivariate(arr_fxy{i});
    
end

% Get degree of polynomials h_{i}(x,y)
vDeg_arr_hxy = vDeg_arr_fxy(1:end-1) - vDeg_arr_fxy(2:end);


% Initialise a cell array
arr_DT1Q1 = cell(nPolys_arr_fxy - 1, 1);

% For each of the polynomials excluding the first f_{1},...,f_{d}
for i = 1 : 1 : nPolys_arr_fxy - 1
    
    
    % Get the degree of f{i}
    m = vDeg_arr_fxy(i+1);
    
    % Get the degree of h_{i}
    n = vDeg_arr_hxy(i);
    
    % Temporarily call the ith entry f(x,y)
    fxy = arr_fxy{i+1};
    
    % Build the matrix T_{n-m}(f(x,y))
    D = BuildD_2Polys(m, n);
    T1 = BuildT1(fxy, m, n);
    Q1 = BuildQ1(n);
    
    arr_DT1Q1{i} = D*T1*Q1;
end

C_fxy = blkdiag(arr_DT1Q1{:});


end



function vRHS = Get_vec_z(arr_zxy)


% Get number of polynomials in array f_{i}(x,y)
nPolys_arr_fxy = size(arr_zxy,1);

% Get the degree of the polynomials f_{i}(x,y)
vDeg_arr_zxy = zeros(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy
    vDeg_arr_zxy(i) = GetDegree_Bivariate(arr_zxy{i});
end


% Initialise an array
arr_rhs = cell(nPolys_arr_fxy,1);

% For each polynomial f_{1},...,f_{d} (Note exclude f_{0})
for i = 1:1:nPolys_arr_fxy
    
    % Get degree of f(x,y)
    m = vDeg_arr_zxy(i);
    
    % Get vector of coefficients of z(x,y)
    v_zxy = GetAsVector(arr_zxy{i});
    
    % Get number of coefficients in z(x,y)
    nCoefficients_zxy = nchoosek(m+2,2);
    
    % Strip zeros from z(x,y)
    v_fxy = v_zxy(1:nCoefficients_zxy);
    
    arr_rhs{i} = v_fxy;
    
end

vRHS = cell2mat(arr_rhs);

end




function C_fxy = BuildDYQ(arr_hxy, vDeg_arr_fxy)
% Build the matrix C(f1,...,fd)
%
% Inputs.
%
% arr_fxy : (Array of Matrices) Array of polynomials f(x,y)
%
% vDeg_arr_fxy : (Vector)
%
%
% Outputs.
%
% C_fxy : (Matrix)

% Get number of polynomials in f_{i}(x,y)
nPolys_arr_hxy = size(arr_hxy,1);

% Get degree of each polynomial f_{i}(x,y)
vDeg_arr_hxy = zeros(nPolys_arr_hxy,1);

for i = 1 : 1 : nPolys_arr_hxy
    
    vDeg_arr_hxy(i) = GetDegree_Bivariate(arr_hxy{i});
    
end


% Initialise a cell array
arr_DT1Q1 = cell(nPolys_arr_hxy,1);

% For each of the polynomials excluding the first f_{1},...,f_{d}
for i = 1:1:nPolys_arr_hxy
    
    
    % Get the degree of h{i}(x,y)
    m = vDeg_arr_hxy(i);
    
    % Get the degree of f{i}(x,y)
    n = vDeg_arr_fxy(i+1);
    
    % Temporarily call the ith entry f(x,y)
    hxy = arr_hxy{i};
    
    % Build the matrix T_{n-m}(f(x,y))
    T1 = BuildT1(hxy,m,n);
    
    D = BuildD_2Polys(m,n);
    Q1 = BuildQ1(n);
    
    arr_DT1Q1{i} = D*T1*Q1;
end

% Include a section of zeros

C_fxy = blkdiag(arr_DT1Q1{:});
nRows = size(C_fxy,1);

nCols = nchoosek(vDeg_arr_fxy(1)+2,2);
zero_section = zeros(nRows,nCols);

C_fxy = [zero_section C_fxy];
end



















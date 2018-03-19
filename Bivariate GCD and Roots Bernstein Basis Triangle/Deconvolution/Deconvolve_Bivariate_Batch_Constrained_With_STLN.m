function [arr_hxy] = Deconvolve_Bivariate_Batch_Constrained_With_STLN(arr_fxy, vDegt_fxy)
%
%
% % Inputs.
%
% arr_fxy : (Array of Matrices) Array of polynomials f(x,y) in Bernstein
% form.
%
% vDegt_fxy : (Vector) Vector of degrees of the polynomials f_{i}(x,y).
%
% % Outputs
%
% arr_hxy : (Array of Matrices)


global SETTINGS

% Get vector of degrees of f_{i}(x,y)
vDegree_arr_fxy = vDegt_fxy;

% Get vector of degrees of h_{i}(x,y) given by f_{i}/f_{i+1}
vDegree_h = diff(vDegree_arr_fxy);

% Get vector of degrees of w_{i}(x,y) given by h_{i}/h_{i+1}
vDegree_wxy = diff([vDegree_h; 0]);

vMultiplicity = find(vDegree_wxy~=0);

% Get the number of polynomials in the array of polynomials f_{i}(x)
nPolys_arr_fxy = size(arr_fxy,1);

% Get the number of polynomials in the array of polynomials h_{i}(x)
nPolys_arr_hxy = nPolys_arr_fxy - 1;

% % Get the degree of each of the polynomials f_{i}(x,y)

% Initialise vector
vDegree_arr_fxy = zeros(nPolys_arr_fxy,1);

% For each polynomial f_{i} get the degree
for i = 1:1:nPolys_arr_fxy
    
    vDegree_arr_fxy(i) = GetDegree_Bivariate(arr_fxy{i});
    
end

% For each polynomial h_{i}
vDeg_arr_hxy = vDegree_arr_fxy(1:end-1) - vDegree_arr_fxy(2:end);


if( SETTINGS.PREPROC_DECONVOLUTIONS)
    
    [th1, th2] = GetOptimalTheta(arr_fxy);
    
else
    th1 = 1;
    th2 = 1;
    
end

% Preprocess the polynomials f_{i}(x,y)
arr_fww = GetPolynomialArrayWithThetas(arr_fxy,th1,th2);


% %
% %
% Build the LHS Matrix
DT_fwwQ = BuildDTQ_2Polys(arr_fww,vDegree_arr_fxy);

% %
% %
% Build the RHS Vector
rhs_fww = BuildRHS_vec(arr_fww);

x_ls = SolveAx_b(DT_fwwQ,rhs_fww);
v_pww = x_ls;
unique_vMult = unique(vMultiplicity);

arr_pww = cell(length(unique_vMult),1);

for i = 1:1:length(unique_vMult)
    
    
    mult = unique_vMult(i);
    
    % Get degree of p(x,y)
    deg_px = vDegt_fxy(mult) - vDegt_fxy(mult + 1);
    
    % Get number of coefficients in p(x,y)
    nCoefficients_px = nchoosek(deg_px + 2, 2);
    
    % Get coefficients of p(x,y) from x_ls
    vec_px = x_ls(1:nCoefficients_px);
    
    % Remove coefficients
    x_ls(1 : nCoefficients_px) =[];
    
    
    nZeros_px = nchoosek(deg_px+1,2);
    
    vec_px = ...
        [
        vec_px;
        zeros(nZeros_px,1)
        ];
    
    arr_pww{i,1} = GetAsMatrix(vec_px,deg_px,deg_px);
    
end

nPolys_arr_pxy = size(arr_pww, 1);
vDeg_arr_pxy = zeros(nPolys_arr_pxy, 1);

for i = 1:1: nPolys_arr_pxy
    
    vDeg_arr_pxy(i,1) = GetDegree_Bivariate(arr_pww{i});
    
end

arr_hww = Get_hxy(arr_pww,unique_vMult);



% %
% %
% Build the array of polynomails z(\omega) which are the structured
% perturbations of the array of polynomials f(x).
arr_zww = cell(nPolys_arr_fxy,1);

for i = 1 : 1 : nPolys_arr_fxy
    
    % Get degree of polynomial f_{i}(x,y)
    m = vDegree_arr_fxy(i);
    
    arr_zww{i} = zeros(m+1,m+1);
    
end

% Build vector q(\omega) consisting of all vectors in z_{i}(x)
v_qww = Get_v_zxy(arr_zww);

% %
% %
% Build the matrix P

% Get number of coefficients across all polynomials f_{i}(x,y)
v_nCoefficients_arr_fxy = zeros(nPolys_arr_fxy,1);

for i = 1:1:nPolys_arr_fxy
    
    m = vDegree_arr_fxy(i);
    v_nCoefficients_arr_fxy(i,1) = nchoosek(m+2,2);
    
end

% Get number of coefficients across all polynomials h_{i}(x,y)
v_nCoefficients_arr_hxy = zeros(nPolys_arr_hxy,1);

for i = 1:1:nPolys_arr_hxy
    
    n = vDeg_arr_hxy(i);
    v_nCoefficients_arr_hxy(i,1) = nchoosek(n+2,2);
    
end

% Get number of coefficients across all polynomials p_{i}(x,y)
v_nCoefficients_arr_pxy = zeros(nPolys_arr_pxy,1);
for i = 1:1:nPolys_arr_pxy
    
    n = vDeg_arr_pxy(i);
    v_nCoefficients_arr_pxy(i,1) = nchoosek(n+2,2);
    
end

% Get number of coefficients in all polynomials f_{i}(x,y)
nCoefficients_fxy = sum(v_nCoefficients_arr_fxy);


nCoefficients_pxy = sum(v_nCoefficients_arr_pxy);

% Get sum of all coefficients except the final f_{i}(x,y)
nCoefficients_rhs = sum(v_nCoefficients_arr_fxy(1:end-1));

% Finally build the matrix P
P = [...
    eye(nCoefficients_rhs) ...
    zeros(nCoefficients_rhs,v_nCoefficients_arr_fxy(end))
    ];


% DY_hQ
DY_hQ = BuildDYQ(arr_hww,vDegree_arr_fxy);

% Set iteration number
ite = 1;

% Build the matrix F
E = [eye(nCoefficients_fxy), zeros(nCoefficients_fxy, nCoefficients_pxy)];

%
%
%
% Build the matrix G

% Build component H_h of G
C_p = DT_fwwQ;
C_q = DY_hQ - P;

C = [C_q C_p];

% Compute the first residual
t = rhs_fww + (P*v_qww) - (DT_fwwQ * v_pww);

% Update Matrix P*z
Pq = P*v_qww;

% Perform test
% for i = 1 : 1 : nPolys_arr_fx
%     vec_fw = [RHS_vec_fww; arr_fx{i}];
% end
%test1 = DY_hQ * vec_fw;
%test2 = DT_fwQ * v_pw;
%test1./test2

condition(ite) = norm(t)./norm(rhs_fww + Pq);

start_point = ...
    [
    v_qww;
    v_pww;
    ];

%Get the iterated value
yy = start_point;

% Get
s = E*(start_point - yy);


while (condition(ite) > SETTINGS.MAX_ERROR_DECONVOLUTIONS)  && ...
        (ite < SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS)
    
    % Use the QR decomposition to solve the LSE problem and then
    % update the solution.
    % min |Fy-s| subject to Gy=t
    y = LSE_new(E, s, C, t);
    
    yy = yy + y;
    
    % output y gives delta p and delta z
    delta_qww = y(1 : nCoefficients_fxy);
    delta_pww = y(nCoefficients_fxy + 1 : nCoefficients_fxy + nCoefficients_pxy);
    
    % Add structured perturbations to vector p(\omega) and z
    v_pww = v_pww + delta_pww;
    v_qww = v_qww + delta_qww;
    
    % Get the updated array of polynomials p_{i}(\omega)
    arr_pww = GetPolynomialArrayFromVector(v_pww, vDeg_arr_pxy);
    arr_zww = GetPolynomialArrayFromVector(v_qww, vDegree_arr_fxy);
    arr_hww = Get_hxy(arr_pww,unique_vMult);
    
    s = E*(start_point - yy);
    
    DY_hQ = BuildDYQ(arr_hww,vDegree_arr_fxy);
    
    % Build the matrix C(f)
    DT_fwwQ = BuildDTQ_2Polys(arr_fww,vDegree_arr_fxy);
    
    % Build the matrix C(z)
    DT_zwwQ = BuildDTQ_2Polys(arr_zww,vDegree_arr_fxy);
    
    % Build component H_h of G
    C_p = DT_fwwQ + DT_zwwQ;
    C_q = DY_hQ - P;

    C = [C_q C_p];
    
    % Update the RHS vector
    rhs_fww = BuildRHS_vec(arr_fww);
    RHS_vec_Pzww = BuildRHS_vec(arr_zww);
    
    
    % Calculate residual and increment t in LSE Problem
    t = ((rhs_fww + RHS_vec_Pzww) - ((DT_fwwQ + DT_zwwQ)*v_pww));
    
    % Increment iteration number
    ite = ite + 1;
    
    % Get condition number
    condition(ite) = norm(t)./norm(rhs_fww + RHS_vec_Pzww);
    
end

if(SETTINGS.PLOT_GRAPHS_DECONVOLUTION)
    
    % Plot termination condiiton
    figure_name = sprintf([mfilename ' : STLN Iterations']);
    figure('name',figure_name)
    plot(log10(condition),'-s','DisplayName','Termination Condition')
    hold off
    
end


% Get array of polynomials h_{i}(x,y)
arr_hxy = GetPolynomialArrayWithoutThetas(arr_hww,th1,th2);



end


function LHS_Matrix = BuildDTQ_2Polys(arr_fxy, vDegt_fxy)
%
% % Inputs
%
% arr_fxy : (Array of Matrices)
%
% % Outputs
%
% LHS_Matrix


vDeg_f = vDegt_fxy;
vDeg_h = diff(vDeg_f);
vDeg_wxy = diff([vDeg_h; 0]);
vMult = find(vDeg_wxy~=0);

% Get number of distinct polynomials h_{i}(x)
nDistinct_hx = length(vMult);

for i = 1:1:nDistinct_hx
    
    if i>1
        old_mult = vMult(i-1);
    else
        old_mult = 0;
    end
    
    new_mult = vMult(i);
    
    arr_Cf{i} = [];
    
    for j = (old_mult+1+1) : 1 : (new_mult+1)
        
        % Get coefficients of previous f(x,y)
        fxy_prev = arr_fxy{j-1};
        
        % Get the degree of previous f(x,y)
        deg_fx_prev = vDeg_f(j-1);
        
        % Get polynomial f(x,y)
        fxy = arr_fxy{j};
        
        % Get the degree of f(x,y)
        deg_fx = vDeg_f(j);
        
        % Get the degree of polynomial h_{i}
        deg_hx = deg_fx_prev - deg_fx;
        
        % Build the cauchy like matrix
        arr_Tf{j} = BuildT1(fxy,deg_fx,deg_hx);
        
        % Build the diagonal matrix D
        arr_D{j} = BuildD_2Polys(deg_fx,deg_hx);
        
        % Stack beneath all other T_f
        arr_Cf{i} = [arr_Cf{i} ; arr_D{j}*arr_Tf{j}];
        
    end
    
    arr_Q{i} = BuildQ1(deg_hx);
    
    DTQ{i} = arr_Cf{i} * arr_Q{i};
end

LHS_Matrix = blkdiag(DTQ{:});

end




function v_zxy =  Get_v_zxy(arr_zxy)
%
% % Inputs
%
%
% arr_zxy : (Array of Matrices)
%
% % Outputs
%
% v_zxy : (Vector)


% Get the number of polynomials in array
nPolys_arr_zxy = size(arr_zxy,1);

% Get the degree of each of the polynomials z(x,y)
vDeg_arr_zxy = zeros(nPolys_arr_zxy);

for i = 1:1:nPolys_arr_zxy
    
    vDeg_arr_zxy(i) = GetDegree_Bivariate(arr_zxy{i});
    
end

% For each polynomial in the array, get as a vector
for i = 1:1:nPolys_arr_zxy
    
    % Get as vector
    v_zxy = GetAsVector(arr_zxy{i});
    
    % Get number of non-zeros
    nZeros_zxy = nchoosek(vDeg_arr_zxy(i)+2,2);
    
    % Remove zeros from vector
    v_zxy = v_zxy(1:nZeros_zxy);
    
    arr_zxy{i} = v_zxy;
end

v_zxy = cell2mat(arr_zxy);

end

function arr_hxy = Get_hxy(arr_pxy, vUniqueMult)
%
% % Inputs
%
% arr_pxy : (Array of Matrices)
%
% vUniqueMult :(Vector)
%
%
% % Outputs
%
% arr_hxy : (Array of Matrices)

% Get number of entries in the array of polynomials p_{i}(x)
nPolys_arr_px = size(arr_pxy,1);

% initialise count
count = 1;


for i = 1:1:nPolys_arr_px
    
    if i == 1
        nReps = vUniqueMult(i);
    else
        nReps = (vUniqueMult(i) - vUniqueMult(i-1));
    end
    
    for j = 1:1:nReps
        arr_hxy{count,1} = arr_pxy{i};
        count = count + 1;
    end
    
end
end


function C_fxy = BuildDYQ(arr_hxy, vDeg_arr_fxy)
% Build the matrix C(f1,...,fd)
%
% Inputs.
%
% arr_hxy : (Array of Matrices) Contains matrices of polynomials f(x,y)
%
% vDeg_arr_fxy : (Vector)
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
arr_DT1Q1 = cell(nPolys_arr_hxy, 1);

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

nColumns = nchoosek(vDeg_arr_fxy(1)+2,2);
zero_section = zeros(nRows,nColumns);

C_fxy = [zero_section C_fxy];
end

function arr_zxy = GetArray(v_zx, v_deg_arr_fxy)
% Given the vector of perturbations of f(x) given by v_zx
%
% % Inputs
%
% v_zx :
%
% v_deg_arr_fxy : (Vector)
%
% % Outputs
%
% arr_zxy :(Array of Matrices)


% Get number of polynomials in arr_fx
nPolys_fxy = size(v_deg_arr_fxy,1);

% Initialise an array
arr_zxy = cell(nPolys_fxy,1);

for i = 1:1:nPolys_fxy
    
    % Get degree of f_{i}(x)
    m = v_deg_arr_fxy(i);
    
    % Get number of coefficients in f_{i}(x)
    nCoefficients_fxy = nchoosek(m+2,2);
    
    % Get the coefficients from the vector
    temp_vec = v_zx(1:nCoefficients_fxy);
    
    % Remove the m+1 coefficients
    v_zx(1:nCoefficients_fxy) = [];
    
    % Get number of zeros in f_{i}(x,y) to be added to form a matrix
    try
        nZeros = nchoosek(m+1,2);
    catch
        nZeros = 0;
    end
    
    % Add zeros
    temp_vec = ...
        [
        temp_vec;
        zeros(nZeros,1)
        ];
    
    % Get as matrix
    mat = GetAsMatrix(temp_vec,m,m);
    
    arr_zxy{i} = mat;
    
    
end


end
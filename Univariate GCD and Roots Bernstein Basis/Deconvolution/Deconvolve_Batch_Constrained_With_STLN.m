function [arr_hx] = Deconvolve_Batch_Constrained_With_STLN(arr_fx, vMult)
% Let vMultiplicities be a vector containing the multiplicities of the
% roots of f_{0}(x). vMult = [m_{1}, m_{2} ,..., m_{n}]
% The division f_{0}/f_{1},...,f_{m_{1}-1} / f_{m_{1}} all have the same solution
%
%
% % Inputs.
%
% arr_fx : (Array of Vectors) Array of polynomials f_{i}(x)
%
% vMult : (Vector) Multiplicities of the factors of f(x) in ascending order
%
% % Outputs.
%
% arr_hx : (Array of Vectors) Array of polynomials h_{i}(x) where
% h_{i}(x) = f_{i-1}(x) / f_{i}(x)

% Global Variables
global SETTINGS

% Get the number of polynomials in the array of polynomials f_{i}(x)
nPolys_arr_fx = size(arr_fx, 1);

% % Get the degree of each of the polynomials f{i}(x)

% Initialise vector
vDeg_arr_fx = zeros(nPolys_arr_fx, 1);

% For each polynomial f_{i} get the degree m_{i}
for i = 1:1:nPolys_arr_fx
    
    vDeg_arr_fx(i) = GetDegree(arr_fx{i});
    
end

% Define M to be the total number of all coefficients of the first d polynomials
% f_{0}...f_{d-1},
M = sum(vDeg_arr_fx + 1) - (vDeg_arr_fx(end) + 1);

% Define M1 to be the total number of all coefficients of polynomials
% f_{0},...,f_{d}
nCoefficients_fx = sum(vDeg_arr_fx + 1);



% Preprocess
if(SETTINGS.PREPROC_DECONVOLUTIONS)
    
    theta = GetOptimalTheta(arr_fx);
    
else
    theta = 1;
    
end

% % Preprocess polynomials f_{i}(x) to get array of polynomials
% f_{i}(\omega)

% Initialise a cell-array for f(w)
arr_fw = GetPolynomialArrayWithThetas(arr_fx, theta);

% % Build LHS Matrix C(f1,...,fd)
DT_fwQ = BuildDTQ(arr_fw, vMult);

% % Build the RHS vector

% RHS vector consists of f_{1},...,f_{m_{i}} where m_{i} is the highest
% degree of any root of f_{0}(x).
RHS_vec_fw = BuildRHSF(arr_fw);


% % Get the vector containing coefficients of polynomials p_{i}(\omega)
v_pw = SolveAx_b(DT_fwQ, RHS_vec_fw);

% Get number of coefficients in all the polynomials p_{i}(\omega)
nCoefficients_px = length(v_pw);

unique_vMult = unique(vMult);

% Get the number of polynomials in the array p_{i}(x)
nPolys_arr_px = length(unique_vMult);

% Get the degree of the polynomials in the array p_{i}(x)
vDegree_arr_px = zeros(nPolys_arr_px, 1);

for i = 1:1:length(unique_vMult)
    
    % Get multiplicity 
    factor_multiplicity = unique_vMult(i);
    
    % Get Degree
    factor_degree = vDeg_arr_fx(factor_multiplicity) - vDeg_arr_fx(factor_multiplicity+1);
    
    % 
    vDegree_arr_px(i) = factor_degree;
end

% %
% %
% Get the polynomials p_{i}(x) repeated to give the set of polynomials
% h_{i}(x).
arr_pw = GetPolynomialArrayFromVector(v_pw, vDegree_arr_px);

% %
% %
% Get the polynomials p_{i}(\omega) repeated to give the set of polynomials
% h_{i}(\omega)
arr_hw = Get_hx(arr_pw, unique_vMult);

% %
% %
% Build the array of polynomails z(\omega) which are the structured
% perturbations of the array of polynomials f(x).
arr_zw = cell(nPolys_arr_fx,1);
for i = 1:1:nPolys_arr_fx
    arr_zw{i} = zeros(vDeg_arr_fx(i) + 1,1);
end

% Build vector z(\omega) consisting of all vectors in z_{i}(x)
v_qw = cell2mat(arr_zw);

% Build the matrix P
P = [eye(M) zeros(M, nCoefficients_fx - M)];

% DY_hQ
DY_hQ = BuildDYQ(arr_hw, vDeg_arr_fx);

% Set iteration number
ite = 1;

% Build the matrix F
E = [eye(nCoefficients_fx) zeros(nCoefficients_fx, nCoefficients_px)];

% Build the matrix G

% Build component H_h of G
% Build component H_z of G
C_z = DY_hQ - P;
C_h = DT_fwQ;

C = [C_z C_h];

% Compute the first residual
t = RHS_vec_fw + (P*v_qw) - (DT_fwQ * v_pw);

% Update Matrix P*z
Pq = P*v_qw;

% Perform test
%for i = 1 : 1 : nPolys_arr_fx
%    vec_fw = [RHS_vec_fw; arr_fx{i}];
%end
%test1 = DY_hQ * vec_fw;
%test2 = DT_fwQ * v_pw;
%test1./test2

vCondition(ite) = norm(t)./norm(RHS_vec_fw + Pq);

start_point = ...
    [
        v_qw;
        v_pw;
    ];

%Get the iterated value
yy = start_point;

% Get
s = E * (start_point - yy);



while (vCondition(ite) > SETTINGS.MAX_ERROR_DECONVOLUTIONS)  && ...
        (ite < SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS)
    
    % Use the QR decomposition to solve the LSE problem and then
    % update the solution.
    % min |Fy-s| subject to Gy=t
    y = LSE_new(E,s,C,t);
    
    yy = yy + y;
    
    % output y gives delta p and delta z
    delta_qw = y(1 : nCoefficients_fx);
    delta_pw = y(nCoefficients_fx + 1 : nCoefficients_fx + nCoefficients_px);
    
    % Add structured perturbations to vector p(\omega)
    v_pw = v_pw + delta_pw;
    v_qw = v_qw + delta_qw;
    
    % Get the updated array of polynomials p_{i}(\omega)
    arr_pw = GetPolynomialArrayFromVector(v_pw, vDegree_arr_px);
    arr_zw = GetPolynomialArrayFromVector(v_qw, vDeg_arr_fx);
    
    arr_hw = Get_hx(arr_pw, unique_vMult);
    
    s = E*(start_point - yy);
    
    DY_hQ = BuildDYQ(arr_hw, vDeg_arr_fx);
    
    % Build the matrix C(f)
    DT_fwQ = BuildDTQ(arr_fw, vMult);
    DT_zwQ = BuildDTQ(arr_zw, vMult);
    
    % Build G
    C_z = DY_hQ - P;
    C_h = DT_fwQ + DT_zwQ;
    
    C = [C_z C_h];
    
    % Update the RHS vector
    RHS_vec_fw = BuildRHSF(arr_fw);
    RHS_vec_Pz = BuildRHSF(arr_zw);
    
    
    % Calculate residual and increment t in LSE Problem
    t = ((RHS_vec_fw+RHS_vec_Pz) - ((DT_fwQ + DT_zwQ)*v_pw));
    
    % Increment iteration number
    ite = ite + 1;
    
    % Get condition number
    vCondition(ite) = norm(t)./norm(RHS_vec_fw + RHS_vec_Pz);
end



% Get array of polynomials h_{i}(x) from h_{i}(\omega)
arr_hx = GetPolynomialArrayWithoutThetas(arr_hw, theta);


%
LineBreakLarge();
fprintf([mfilename ' : ' sprintf('Required Number of iterations : %i \n',ite)]);
LineBreakLarge();

if (SETTINGS.PLOT_GRAPHS_DECONVOLUTION_LRA)
    figure()
    hold on
    plot(log10(vCondition))
    hold off
end

end

function LHS_Matrix = BuildDTQ(arr_fx, vMult)
% %
% %
% Build the LHS Coefficient matrix

% For each distinct hx, build the partition of the Matrix C(f_{m_{i}+1},...,f_{m_{i+1}})
% C(f_{1},...f_{m1}), C(f_{m1+1},...,f_{m2}),...
% C(f_{1},...,f_{m1}) = [T(f1) ; T(f2) ; ... T(fm1)]

% Get number of distinct polynomials in h_{i}
nDistinct_hx = length(vMult);

arr_Cf = cell(nDistinct_hx,1);
arr_Q = cell(nDistinct_hx,1);
arr_DTQ = cell(nDistinct_hx,1);
for i = 1:1:nDistinct_hx
    
    if i > 1
        old_mult = vMult(i-1);
    else % No previous multiplicity
        old_mult = 0;
    end
    
    new_mult = vMult(i);
    
    arr_Cf{i} = [];
    
    % for each polynomial f_{i} in the interval f_{m_{i-1}+1}...f_{m_{i}}
    
    arr_Tf = cell((new_mult+1) - (old_mult+1+1));
    arr_D = cell((new_mult+1) - (old_mult+1+1));
    
    for j = (old_mult+1+1) : 1 : (new_mult+1)
        
        % Get the degree of the previous polynomial f_{i-1}(x)
        fx_prev = arr_fx{j-1};
        deg_fx_prev = GetDegree(fx_prev);
        
        % Get the degree of the current polynomial f_{i}(x)
        fx = arr_fx{j};
        deg_fx = GetDegree(fx);
        
        % Get the degree of the polynomial h_{i} = f_{i-1}/f_{i}
        deg_hx = deg_fx_prev - deg_fx;
        
        % Build the Cauchy like matrix T_{m_{i} - m_{i-1}}(f_{i})
        arr_Tf{j,1} = BuildT1(fx, deg_hx);
        
        arr_D{j,1} = BuildD_2Polys(deg_fx, deg_hx);
        
        % Stack beneath all other T_{f} which are multiplied by [_{i}(x)
        arr_Cf{i,1} = [arr_Cf{i} ; arr_D{j}*arr_Tf{j}];
    end
    
    arr_Q{i} = BuildQ1(deg_hx);
    
    arr_DTQ{i} = arr_Cf{i} * arr_Q{i};
    
end


LHS_Matrix = blkdiag(arr_DTQ{:});

end




function arr_hx = Get_hx(arr_px, vUniqueMult)
%
% % Inputs
%
% arr_px : (Array of Vectors) Array of vectors containing coefficietns of 
% the polynomials p_{i}(x)
% 
% vUniqueMult : (Vector)
%
% % Outputs
%
% arr_hx : (Array of Vectors) Array of polynomials h_{i}(x)

% Get number of entries in the array of polynomials p_{i}(x)
nEntries_arr_px = size(arr_px,1);

% initialise count
count = 1;


for i = 1 : 1 : nEntries_arr_px
    
    if i == 1
        nRepititions = vUniqueMult(i);
    else
        nRepititions = (vUniqueMult(i) - vUniqueMult(i-1));
    end
    
    % Insert the vector n times
    for j = 1 : 1 : nRepititions
        arr_hx{count, 1} = arr_px{i};
        count = count + 1;
    end
    
end
end





function Y_new = BuildDYQ(arr_hw,vDeg_arr_fx)
% Build the coefficient matrix DYU. This is the change of variable such
% that
% D^{-1}*E(z)*Q * g = D^{-1}*Y(g)*U * z

% Inputs.
%
% arr_hw : Set of polynomials h_{i}(w)
%
% vDeg_arr_fx : vector of degrees of polynomials f_{0},...

nPolys_hw = size(arr_hw,1);

for i = 1:1:nPolys_hw
    
    % Start with f1*h1
    % h_{1} is the first in the cell array h_{i}
    % f_{1} is the second in the cell array f_{i}
    % deg(f_{1}) = m(2)
    
    % Get polynomial h(w)
    hw = arr_hw{i,1};
    
    % Get degree of f_{i}
    deg_fw = vDeg_arr_fx(i+1);
    
    y{i,1} = real(BuildD0Y1U1(hw,deg_fw));
end

%Build the Coefficient Matrix C
num_Rows = 0;
for i = 1:length(vDeg_arr_fx)-1
    num_Rows = num_Rows + 1 + (vDeg_arr_fx(i));
end
cols = (vDeg_arr_fx(1)+1);

xx = zeros(num_Rows,cols);
Y = blkdiag( y{1:length(y)});
Y = [xx Y];

Y_new = Y;


end

function Y1 = BuildD0Y1U1(hx,m1)
global SETTINGS

if( SETTINGS.BOOL_LOG)
    
    Y1 = BuildD0Y1U1_log(hx,m1);
    
else
    Y1 = BuildD0Y1U1_nchoosek(hx,m1);
    
end
end

function Y1 = BuildD0Y1U1_nchoosek(hx,m1)
% Build the Partition of the Coefficient matrix D_{i-1}Y_{i}U_{i}
% to perform the multiplication C(h_{i}) * f_{i} = f_{i+1}
%
% Inputs.
%
% hx : Coefficients of polynomial h_{i}(x)
%
% m1 : Degree of polynomial f_{i}



% Get degree of polynomial h(x) where deg(h_{1}) = n_{1} = m_{0} - m_{1}
n1 = GetDegree(hx);

% Y1 = zeros(m0+1,m1+1);
Y1 = [];

% for each column j = 1:1:m0-m1+1
for j = 0:1:m1
    % for each row i = 1:1:m1+1
    for i = j:1:j+n1
        Y1(i+1,j+1) = ...
            hx(i-j+1) .* nchoosek(i,j) .* nchoosek(n1+m1-i,m1-j);
        
    end
end


Y1 = Y1 ./  nchoosek(n1+m1,m1);

end


function Y1 = BuildD0Y1U1_log(hx,m1)
% Build the Partition of the Coefficient matrix D_{i-1}Y_{i}U_{i}
% h
% m0 = degree of previous polynomial
% m1 = degree of current polynomial


% Get Degree of polynomial deg(h_{x}) = n_{1} = m_{0}-m_{1}
n1 = GetDegree(hx);

m0 = n1+m1;


Y1 = [];
% for each column i = 1:1:m0-m1+1
for k = 0:1:m1
    % for each row j = 1:1:m1+1
    for j = k:1:k+n1
        BinomsEval_Log = lnnchoosek(j,k) + lnnchoosek(m0-j,m1-(j-k));
        BinomsEval_Exp = 10.^BinomsEval_Log;
        Y1(j+1,k+1) = hx(j-k+1) .* BinomsEval_Exp;
    end
end

% Include the denominator
Denom_Log = lnnchoosek(n1+m1,m1);
Denom_Exp = 10.^Denom_Log;

Y1 = Y1 ./  Denom_Exp;


end

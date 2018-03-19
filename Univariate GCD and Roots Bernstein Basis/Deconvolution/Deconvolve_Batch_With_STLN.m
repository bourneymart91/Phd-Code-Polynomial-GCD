function [arr_hx] = Deconvolve_Batch_With_STLN(arr_fx)
% Performs a series of d deconvolutions over a set of polynomials,
% where each polynomial g_{i} appears in two of the deconvolutions.
%
%
% Input
%
% arr_fx : (Array of Vectors) Each cell of the array contains the coefficients 
% of the polynomials f_{i}(x) 
%
% Output:
%
% arr_hx : (Array of Vectors) Vectors containing coefficients of polynomials
% h_{i}(x) where h_{i} = f_{i}/f_{i+1}



% Global Variables
global SETTINGS

% Get the number of polynomials in the array of f_{i}(x)
nPolynomials_arr_fx = size(arr_fx, 1);

% % Get the degree m_{i} of each of the polynomials f_{i}

% Initialise vector to store degrees of f_{i}(x)
vDegree_arr_fx = zeros(nPolynomials_arr_fx, 1);

% For each polynomial f_{i}, get its degree.
for i = 1 : 1 : nPolynomials_arr_fx
    
    vDegree_arr_fx(i) = GetDegree( arr_fx{i});
    
end

% % Get the degrees n{i} of polynomials h_{i} = f_{i-1}/f_{i}.
vDegree_arr_hx = (vDegree_arr_fx(1:end-1) - vDegree_arr_fx(2:end));


% Define M to be the total number of coefficeints of all the polynomials
% f_{i} excluding the last f_{i}.
% f_{0},...,f_{nPolys_f-1}.
M = sum(vDegree_arr_fx + 1) - (vDegree_arr_fx(end:end)+1);

% Define M1 to be the total number of all coefficients of all of the polynomials
% f_{i}
nCoefficients_fx = sum(vDegree_arr_fx + 1);

% Define N to be the number of coefficients of all h_{i}
nCoefficients_hx = sum(vDegree_arr_hx + 1);

% Preprocess the polynomials to get optimal value \theta
if(SETTINGS.PREPROC_DECONVOLUTIONS)
    
    theta = GetOptimalTheta(arr_fx);
    
else
    
    theta = 1;
    
end

% Preprocess polynomials f_{i}(x)
arr_fw = GetPolynomialArrayWithThetas(arr_fx, theta);

% % Write Deconvolutions in form [D^{-1}T(f)Q] h = RHS_f

% Get the right hand side vector of coefficients of f_{\omega}
vRHS_fw = BuildRHSF(arr_fw);

% Get the Left hand side matrix C(f1,...,fd)
DTQ = BuildDTQ(arr_fw);

% Get the solution vector h(w) in the system of equations
% DCQ * hw = RHS_vec.
v_hw = SolveAx_b(DTQ, vRHS_fw);


% Separate solution vector h, into component parts h_{1},h_{2},...h_{d},
% each of degree n_{i}

% GetArray of polynomials h_{i}(\omega)
arr_hw = GetPolynomialArrayFromVector(v_hw, vDegree_arr_hx);

% Let z be  vectors of perturbations to polynomials f_{i} such that
% z = [z{0} z{1} z{2} z{3} ... z{d}]
arr_qw = cell(nPolynomials_arr_fx, 1);

for i = 1 : 1 : nPolynomials_arr_fx
    
    % initialise polynomial z_{i} as a zero vector.
    arr_qw{i,1} = zeros(vDegree_arr_fx(i) + 1,1);
    
    
end

% Build vector z, consisting of all vectors z_{i}
v_qw = cell2mat(arr_qw);

% Build the Matrix P
P = [eye(M) zeros(M, nCoefficients_fx - M)];

% Get Vector of perturbations for RHS by multiplying perturbation vector by
% P, such that we eliminate the z_max

% Build Matrix Y, where E(z)h = Y(h)z
DYU = BuildDYU(arr_hw, vDegree_arr_fx);

% Compute the initial residual
t = (vRHS_fw + (P*v_qw) - (DTQ*v_hw));

% Set the iteration counter.
ite = 1;

% Initialise the weight vector
E = [eye(nCoefficients_fx) zeros(nCoefficients_fx, nCoefficients_hx)];

C_q = (DYU) - P;
C_h = DTQ;
C = [C_q C_h];



condition(ite) = norm(t)./ norm(vRHS_fw);

start_point = ...
    [
        v_qw;
        v_hw;
    ];

yy = start_point;

s = E*(start_point - yy);
% Perform iteration to obtain perturbations

while (condition(ite) > SETTINGS.MAX_ERROR_DECONVOLUTIONS)  && ...
        (ite < SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS)
    
    % Use the QR decomposition to solve the LSE problem and then
    % update the solution.
    % min |Fy-s| subject to Gy=t
    y = LSE(E, s, C, t);
    
    % Update vector yy
    yy = yy + y;
    
    % Output y gives delta h and delta z
    delta_z = y(1 : nCoefficients_fx );
    delta_h = y( nCoefficients_fx + 1 : nCoefficients_fx + nCoefficients_hx);
    
    % Add structured perturbations to vector h(\omega).
    v_qw = v_qw + delta_z;
    v_hw = v_hw + delta_h;
    
    
    % Separate delta_z into its component vectors delta_z0 delta_z1,...,
    % delta_zd
    arr_qw = GetPolynomialArrayFromVector(v_qw, vDegree_arr_fx);
    arr_hw = GetPolynomialArrayFromVector(v_hw, vDegree_arr_hx);
    
    
    %Increment s in LSE Problem
    s = E*(start_point - yy);
    
    %Build iterative DYU
    DYU = BuildDYU(arr_hw, vDegree_arr_fx);
    
    % Build DCEQ
    DC_fQ = BuildDTQ(arr_fw);
    DC_zQ = BuildDTQ(arr_qw);
    
    C_q = (DYU) - P;
    C_h = DC_fQ + DC_zQ;
    C = [C_q C_h];
    
    % Update the RHS_vector
    vRHS_fw = BuildRHSF(arr_fw);
    vRHS_zw = BuildRHSF(arr_qw);
    
    % Calculate residual and increment t in LSE Problem
    t = ((vRHS_fw + vRHS_zw ) - ((DC_fQ + DC_zQ)*v_hw));
    
    % Get the condition
    condition(ite +1) = norm(t)./norm((vRHS_fw + vRHS_zw));
    
    % Increment iteration number
    ite = ite + 1;
    
end % End of loop

% Get array of polynomials h_{i}(x) from h_{i}(\omega)
arr_hx = GetPolynomialArrayWithoutThetas(arr_hw, theta);



% Print outputs to command line
LineBreakLarge();
fprintf([mfilename ' : ' sprintf('Required Number of iterations : %i \n',ite)]);
LineBreakLarge();

if (SETTINGS.PLOT_GRAPHS_DECONVOLUTION_LRA)
    fig_name = sprintf([mfilename ' : ' 'Condition at Iterations']);
    figure('name',fig_name)
    hold on
    plot(log10(condition),'-s','DisplayName','Condition')
    hold off
end

end

function Y_new = BuildDYU(arr_hx, vDeg_arr_fx)
% Build the coefficient matrix DYU. This is the change of variable such
% that
% D^{-1}*E(z)*Q * g = D^{-1}*Y(g)*U * z

% Inputs.
%
% arr_hw : (Array of Vectors) Set of polynomials h_{i}(x)
%
% vDeg_arr_fx : (Vector) Vector of degrees of polynomials f_{0},...

nPolys_arr_hx = size(arr_hx,1);

for i = 1:1:nPolys_arr_hx
    
    % Start with f1*h1
    % h_{1} is the first in the cell array h_{i}
    % f_{1} is the second in the cell array f_{i}
    % deg(f_{1}) = m(2)
    
    % Get polynomial h(w)
    hx = arr_hx{i,1};
    
    % Get degree of f_{i}
    deg_fx = vDeg_arr_fx(i+1);
    
    y{i,1} = real(BuildD0Y1U1(hx, deg_fx));
end

%Build the Coefficient Matrix C
nRows = 0;
for i = 1:length(vDeg_arr_fx)-1
    nRows = nRows + 1 + (vDeg_arr_fx(i));
end
nCols = (vDeg_arr_fx(1)+1);

xx = zeros(nRows,nCols);
Y = blkdiag( y{1:length(y)});
Y = [xx Y];

Y_new = Y;


end

function Y1 = BuildD0Y1U1(hx, m1)
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
% hx : (Vector) Coefficients of polynomial h_{i}(x)
%
% m1 : (Int) Degree of polynomial f_{i}

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


function Y1 = BuildD0Y1U1_log(hx, m1)
% Build the Partition of the Coefficient matrix D_{i-1}Y_{i}U_{i}
% h
%
% % Inputs
%
% hx : (Vector) Vector of Coefficients of h(x)
%
% m1 : (Int) degree of current polynomial
%
% % Outputs
%
% Y1 : (Matrix)

% Get Degree of polynomial deg(h_{x}) = n_{1} = m_{0}-m_{1}
n1 = GetDegree(hx);

m0 = n1+m1;


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









function DCQ = BuildDTQ(arr_fx)
% Build the left hand side convolution matrix for the set of deconvolutions
%
%
% Inputs
%
% arr_fx : (Array of Vectors) Array of vectors containing coefficients of
% the polynomials f_{i}(x)

% Get number of polynomials in array of f_{i}(x)
nPolys_arr_fx = size(arr_fx,1);

% For each of the polynomials f_{i}(x), excluding the final polynomial
arr_DT1Q1 = cell(nPolys_arr_fx -1);

for i = 2:1:nPolys_arr_fx
    
    % Get the polynomial f_{i} = set_f{i+1}
    fx = arr_fx{i};
    
    % Get the polynomial f_{i-1} = set_f{i}
    fx_prev = arr_fx{i-1};
    
    % Get degree of polynomial f_{i} = m_{i}
    deg_fx = GetDegree(fx);
    
    % Get the degree of polynomial f_{i-1}
    deg_fx_prev = GetDegree(fx_prev);
    
    % Get the degree of the polynomial h_{i}
    deg_hw = deg_fx_prev - deg_fx;
    
    % Build the Matrix T(f)
    T1 = BuildT1(fx, deg_hw);
    
    % Build the matrix D^{-1}_{}
    D = BuildD_2Polys(deg_fx, deg_hw);
    
    % Build the matrix Q_{}
    Q1 = BuildQ1(deg_hw);
    
    % Add the matrix DTQ to an array of matrices ready to construct the LHS
    % matrix
    arr_DT1Q1{i-1}  = D*T1*Q1;
    
end

%Build the Coefficient Matrix C of all matrices c
DCQ = blkdiag(arr_DT1Q1{1:length(arr_DT1Q1)});

end


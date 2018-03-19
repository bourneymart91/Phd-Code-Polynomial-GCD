function arr_hx = Deconvolve_Batch_With_STLN(arr_fx)
% Given the set of polynomials f_{0},...,f_{1}. compute the series of
% deconvolutions h_{1} = f_{1}/f_{0} h_{2} = f_{2}/f_{1},...
% Perform the deconvolutions by producing the structure
% diag [ C(f_{1}) C(f_{2}) ... ] [h1 h2 ...]^{T} = [f_{0} f_{1} ...]
%
% % Inputs
%
% arr_fx : (Array of Vectors) Array of vectors containing coefficients of
% the polynomials f_{i}(x)
%
%
% % Outputs
%
% arr_hx : (Array of Matrices) Array of vectors containing coefficients of
% the polynomials h_{i}(x)

% Global Variables.
global SETTINGS


% Get the number of polynomials in the set set_f
nPolys_fx = size(arr_fx, 1);

% Get the degree m_{i} of each of the polynomials f_{i} and store in a
% vector.

vDeg_arr_fx = zeros(nPolys_fx, 1);
for i = 1:1:nPolys_fx
    vDeg_arr_fx(i) = GetDegree(arr_fx{i});
end

% Get the degrees n{i} of polynomials h_{i} = f_{i}/f_{i+1}.
vDeg_arr_hx = (vDeg_arr_fx(1:end-1) - vDeg_arr_fx(1+1:end))';


% Define M to be the total number of all coefficients of the first d polynomials
% f_{0}...f_{d-1},
M = sum(vDeg_arr_fx+1) - (vDeg_arr_fx(end)+1);

% Define M1 to be the total number of all coefficients of polynomials
% f_{0},...,f_{d}
nCoefficients_fx = sum(vDeg_arr_fx+1);

% Define N to be the number of coefficients of all h_{i}
nCoefficients_hx = sum(vDeg_arr_hx+1);


% % Preprocess

if(SETTINGS.PREPROC_DECONVOLUTIONS)
    
    theta = GetOptimalTheta(arr_fx,vDeg_arr_fx);
    
else
    theta = 1;
end

% Get polynomials f_{i}(\omega) from f_{i}(x) 
arr_fw = GetPolynomialArrayWithThetas(arr_fx, theta);


% Build the LHS Matrix C(f1,...fd)
Cf = BuildC(arr_fw);

% Build the RHS vector

% RHS vector consists of f_{1},...,f_{m_{i}} where m_{i} is the highest
% degree of any root of f_{0}(x).
RHS_vec_f = BuildRHSF(arr_fw);

% Solve h_{0} for initial values of h
v_hw = SolveAx_b(Cf, RHS_vec_f);

% Build the array of polynomials h_{i}(x)
arr_hw = GetPolyArrayFromVector(v_hw, vDeg_arr_hx);

% Let z be  vectors of perturbations to polynomials fi such that
% z = [z{0} z{1} z{2} z{3} ... z{d}]
arr_zw = cell(1,nPolys_fx);

for i = 1 : 1 : nPolys_fx
    arr_zw{i} = zeros(1,vDeg_arr_fx(i)+1);
end

% Build vector z, consisting of all vectors z_{i}
z_o = [arr_zw{:}]';

% Build the Matrix P
P = [eye(M) zeros(M,nCoefficients_fx - M)];

% Build Matrix Y, where E(z)h = Y(h)z
Y_h = BuildY(arr_hw,vDeg_arr_fx);

% Set the iteration counter.
ite = 1;

% Build the identity matrix F.
F = eye(nCoefficients_hx + nCoefficients_fx);

% % Build the matrix G

H_h = Cf;

H_z = Y_h - P;

G = [H_h H_z];

% %
% %
% Compute the first residual
res_vec = (RHS_vec_f + (P*z_o) - (Cf*v_hw));

% Update Matrix Pz
v_zw = z_o;

Pz = P*v_zw;

% Get the initial residual
condition(ite) = norm(res_vec)./norm(RHS_vec_f + Pz);

% Get the start point aka : y^(0)
start_point = ...
    [
    v_hw;
    z_o;
    ];

% Get the iterated value = y^(j)
yy = start_point;

% Get -y^(j)-y^(0)
s = -(yy - start_point);

% %
% %
% Perform iteration to obtain perturbations

while (condition(ite) > SETTINGS.MAX_ERROR_DECONVOLUTIONS)  && ...
        (ite < SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS)
    
    % Use the QR decomposition to solve the LSE problem and then
    % update the solution.
    % min |Fy-s| subject to Gy=t
    y = LSE(F, s, G, res_vec);
    
    % Let yy be the vector of latest version of (h+\delta h), (z + \delta z)
    yy = yy + y;
    
    % Output y gives delta h and delta z
    delta_h = y(1 : nCoefficients_hx);
    delta_z = y(nCoefficients_hx+1 : end);
    
    % Add structured perturbations to vector hx.
    v_hw = v_hw + delta_h;
    
    % Add structured perturbations to vector z.
    v_zw = v_zw + delta_z;
    
    % Get the updated array of polynomials h_{i}(x)
    arr_hw = GetPolyArrayFromVector(v_hw, vDeg_arr_hx);
    
    % Separate zx into an array
    arr_zw = GetPolyArrayFromVector(v_zw, vDeg_arr_fx);
    
    % Increment s in LSE Problem
    s = -(yy - start_point);
    
    % Build iterative Y, where Y
    Y_h = BuildY(arr_hw, vDeg_arr_fx);
    
    % Build the matrix C(f,...,f)
    Cf = BuildC(arr_fw);
    
    % Build the matrix C(z,...,z)
    Cz = BuildC(arr_zw);
    
    % % Build Matrix G for LSE problem
    
    H_z = Y_h - P;
    
    H_h = Cf + Cz;
    
    G = [H_h H_z];
    
    % Update the RHS Vector
    RHS_vec_f   = BuildRHSF(arr_fw);
    RHS_vec_Pz  = BuildRHSF(arr_zw);
    
    % Calculate residual and increment t in LSE Problem
    res_vec = (RHS_vec_f + RHS_vec_Pz) - ((Cf+Cz)*v_hw);
    
    % Increment iteration number
    ite = ite + 1;
    
    % Get condition number
    condition(ite) = norm(res_vec) ./ norm(RHS_vec_f + RHS_vec_Pz);
    
    
end



% Print outputs to command line
fprintf([mfilename ' : ' 'Performed Deconvolutions\n'])
fprintf([mfilename ' : ' sprintf('Iterations required for Batch Deconvolution %i\n', ite)])

if(SETTINGS.PLOT_GRAPHS)
    
    figure_name = sprintf('%s : Condition',mfilename);
    figure('name',figure_name)
    hold on
    plot(log10(condition),'-s')
    hold off
    
end

% Get polynomials h_{i}(x) from h_{i}(\omega)
arr_hx = GetPolynomialArrayWithoutThetas(arr_hw,theta);


end




function C = BuildC(arr_fx)
% Build the matrix C for the series of polynomial deconvolution
% f1/f2 = h1, f2/f3 = h2, ... , where the solutions h_{i} are given in the
% vector x of the expression Cx = b.
% The matrix  is a strucutred matrix of coefficients of polynomials
% f_{2},...,f_{n} where n is the number of polynomials in the set f_{i}. C
% has a block diagonal structure, where each block T_{i} on the diagonal is a
% toeplitz matrix of coefficients of f_{i+1}
%
% % Inputs.
%
% arr_fx = (Array of Vectors) array of polynomials f(x)
%
% % Outputs.
%
% C : (Matrix) Matrix DCQ containing coefficients of f_{i} as described above.

% Get the number of polynomials in the array f_{i}
nPolys_arr_fx = length(arr_fx);

% Initialise an array of cells to store the submatrices which will make up
% DCQ
arr_T = cell(nPolys_arr_fx - 1, 1);

% For each polynomial f_{i}
for i = 1 : 1 : (nPolys_arr_fx-1)
    
    % Get the polynomial f_{i} = set_f{i+1}
    fw = arr_fx{i+1};
    
    % Get the polynomial f_{i-1} = set_f{i}
    fw_prev = arr_fx{i};
    
    % Get degree of polynomial f_{i} = m_{i}
    deg_fw = GetDegree(fw);
    
    % Get the degree of polynomial f_{i-1}
    deg_fw_prev = GetDegree(fw_prev);
    
    % Get the degree of the polynomial h_{i}
    deg_hw = deg_fw_prev - deg_fw;
    
    % Build the Matrix T(f)
    arr_T{i} = BuildT1(fw,deg_hw);
end

% Build the Coefficient Matrix C of all matrices c
C = blkdiag(arr_T{1:length(arr_T)});

end

function Y = BuildY(arr_hx, vec_m)
% Build the coefficient matrix Y. This is the change of variable such
% that
% E(z)*h = Y(h)*f
%
% % Inputs
%
% arr_hx : (Array of Vectors) Contains coefficients of polynomials
% h_{i}(x)
%
% vec_m : (Vector) Degree of each polynomial h_{i}(x)
%
% % Outputs
%
% Y : (Matrix)

% Get number of polynomials in array of h_{i}(x)
nPolys_arr_hx = size(arr_hx,1);

% Initialise cell array of matrices C(h(x))
arr_Ch = cell(nPolys_arr_hx,1);

for i = 1:1:nPolys_arr_hx
    
    % Start with f1*h1
    % h_{1} is the first in the cell array h_{i}
    % f_{1} is the second in the cell array f_{i}
    % deg(f_{1}) = m(2)
    
    % Get polynomial h_{i}(x)
    hx = arr_hx{i};
    
    % Get degree of f_{i}(x)
    deg_fx = vec_m(i+1);
    
    % build the matrix C_{m_{i}}(h_{i}(x))
    arr_Ch{i} = (BuildY1(hx,deg_fx));
end

% Build the Coefficient Matrix C

% Must include a diagonal section of zeros, since C*f =
nRows = 0;
for i = 1:length(vec_m)-1
    nRows = nRows + 1 + (vec_m(i));
end
cols = (vec_m(1)+1);

xx = zeros(nRows,cols);
Y = blkdiag( arr_Ch{1:length(arr_Ch)});
Y = [xx Y];



end


function Y1 = BuildY1(hx,m1)
% Construct a partition Y1 of the matrix Y.
Y1 = BuildT1(hx,m1);
end



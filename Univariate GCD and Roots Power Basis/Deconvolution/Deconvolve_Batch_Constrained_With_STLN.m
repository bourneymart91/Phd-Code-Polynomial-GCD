function [arr_hx] = Deconvolve_Batch_Constrained_With_STLN(arr_fx, vMult)
% Let vMultiplicities be a vector containing the multiplicities of the
% roots of f_{0}(x). vMult = [m_{1}, m_{2} ,..., m_{n}]
% The division f_{0}/f_{1},...,f_{m_{1}-1} / f_{m_{1}} all have the same solution
%
%
% % Inputs.
%
% arr_fx : *(Array of Vectors) Array containing coefficients of polynomials
% f_{i}(x)
%
% vMult : (Vector) Multiplicities of the factors of f(x) in ascending order.
%
% % Outputs.
%
% arr_hx : Array of polynomials h(x) where h_{i}(x) = f_{i-1}(x) / f_{i}(x)


% Global Variables.
global SETTINGS


% Get the number of polynomials in the array of f_{i}(x)
nPolys_arr_fx = length(arr_fx);

% % Get the degree of each polynomial 

% Initialise vector to store degree of polynomials f_{i}(x)
vDeg_arr_fx = zeros(nPolys_arr_fx,1);

% For each polynomial f_{i} get the degree
for i = 1:1:nPolys_arr_fx
    vDeg_arr_fx(i) = GetDegree(arr_fx{i});
end

% Define M to be the total number of all coefficients of the first d polynomials
% f_{0},...,f_{d-1},
M = sum(vDeg_arr_fx+1) - (vDeg_arr_fx(end)+1);

% Define M1 to be the total number of all coefficients of polynomials
% f_{0},...,f_{d}
nCoefficients_fx = sum(vDeg_arr_fx + 1);

% % Preprocess polynomials f_{i}(x) to obtain f_{i}(\omega)

% Get optimal value of theta
if(SETTINGS.PREPROC_DECONVOLUTIONS)
    theta = GetOptimalTheta(arr_fx, vDeg_arr_fx);
else
    theta = 1;
end

arr_fw = GetPolynomialArrayWithThetas(arr_fx, theta);


% % Build LHS Matrix C(f1,...,fd)
C_fw = BuildC(arr_fw, vMult);

% %
% %
% Build the RHS vector

% RHS vector consists of f_{1},...,f_{m_{i}} where m_{i} is the highest
% degree of any root of f_{0}(x).
RHS_vec_fw = BuildRHSF(arr_fw);

% Get vector of coefficients of h_{i}(x) for all i.
v_pw = SolveAx_b(C_fw,RHS_vec_fw);

nCoefficients_px = length(v_pw);

% Get unique multiplities from the multiplicity vector
unique_vMult = unique(vMult);

nPolys_arr_px = length(unique_vMult);


vDeg_arr_px = zeros(nPolys_arr_px , 1);

for i = 1:1: nPolys_arr_px
    mult = unique_vMult(i);
    deg = GetDegree(arr_fx{mult}) - GetDegree(arr_fx{mult+1});
    vDeg_arr_px(i) = deg;
end

% %
% %
% Get the polynomials p_{i}(x) repeated to give the set of polynomials
% h_{i}(x).
arr_pw = GetPolyArrayFromVector(v_pw,vDeg_arr_px);

% %
% %
% Get the polynomials p_{i}(x) repeated to give the set of polynomials
% h_{i}(x).
arr_hw = Get_hx(arr_pw,unique_vMult);



% Build the array of polynomials z(x) which are the structured
% perturbations of the array of polynomials f(x).


arr_zw = cell(nPolys_arr_fx, 1);
for i = 1 : 1 : nPolys_arr_fx
    arr_zw{i} = zeros(vDeg_arr_fx(i) +1,1);
end

% Build the vector zx consisting of all vectors in arr_zx
v_zw = cell2mat(arr_zw);

% Build the matrix P
P = [eye(M) zeros(M,nCoefficients_fx-M)];

% Build Matrix Y, where E(z)h = Y(h)z
Y_h = BuildY(arr_hw, vDeg_arr_fx);

% Set the iteration number
ite = 1;

% Build the identity matrix F.
F = eye(nCoefficients_px + nCoefficients_fx);

% %
% %
% Build the matrix G

% Build component H_h of G
H_h = C_fw;

% Build component H_z of G
H_z = Y_h - P;

G = [H_h H_z];

% %
% %
% Compute the first residual
res_vec = RHS_vec_fw + (P*v_zw) - (C_fw * v_pw);

% Update Matrix P*z
Pz = P*v_zw;



% Perform test
%--------------
%for i = 1 : 1 : nPolys_arr_fx
%    vec_fw = [RHS_vec_fw; arr_fx{i}];
%end
%test1 = Y_h*vec_fw;
%test2 = C_fw*v_pw;
%test1./test2
%--------------


% Get the intial residual
condition(ite) = norm(res_vec)./norm(RHS_vec_fw + Pz);

% Get the start point aka : y^(0)
start_point = ...
    [
    v_pw;
    v_zw;
    ];

% Get the iterated value = y^(j)
yy = start_point;

% Get -y^(j)-y^(0)
s = -(yy - start_point);

% %
% %
% Perform iteration to obtain perturbations

while (condition(ite) > SETTINGS.MAX_ERROR_DECONVOLUTIONS) && ...
        (ite < SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS)
    
    % Use the QR decomposition to solve the LSE problem and then
    % update the solution.
    % min |Fy-s| subject to Gy=t
    y = LSE(F,s,G,res_vec);
    
    yy = yy + y;
    
    % output y gives delta p and delta z
    delta_pw = y(1:nCoefficients_px);
    delta_zw = y(nCoefficients_px+1:end);
    
    % Add structured perturbations to vector hx.
    v_pw = v_pw + delta_pw;
    
    % Add structured perturbations to vector z(\omega).
    v_zw = v_zw + delta_zw;
    
    % Get the updatated array of polynomials p_{i}(x)
    arr_pw = GetPolyArrayFromVector(v_pw,vDeg_arr_px);
    
    % Get the updated array of polynomials of z_{i}(x)
    arr_zw = GetPolyArrayFromVector(v_zw,vDeg_arr_fx);
    
    % Get the updated array of polynomials h_{i}(x)
    arr_hw = Get_hx(arr_pw,unique_vMult);
    
    % Increment s in LSE Problem
    s = -(yy - start_point);
    
    % Build the matrix Y(h_{i})
    Y_h = BuildY(arr_hw,vDeg_arr_fx);
    
    % Build the matrix C(f)
    C_fw = BuildC(arr_fw,vMult);
    
    % Build the matrix C(z)
    C_zw = BuildC(arr_zw,vMult);
    
    % Build G
    H_z = Y_h - P;
    H_h = C_fw + C_zw;
    
    G = [H_h H_z];
    
    % Update the RHS vector
    RHS_vec_fw = BuildRHSF(arr_fw);
    RHS_vec_Pz = BuildRHSF(arr_zw);
    
    
    % Calculate residual and increment t in LSE Problem
    res_vec = ((RHS_vec_fw + RHS_vec_Pz) - ((C_fw + C_zw)*v_pw));
    
    % Increment iteration number
    ite = ite + 1;
    
    % Get condition number
    condition(ite) = norm(res_vec)./norm(RHS_vec_fw + RHS_vec_Pz);
    
    
end

% Get h_{i}(x) from h_{i}(\omega)
arr_hx = GetPolynomialArrayWithoutThetas(arr_hw, theta);

% Print outputs to command line
LineBreakLarge()
fprintf([mfilename ' : ' 'Performed Deconvolutions \n'])
fprintf([mfilename ' : ' sprintf('Iterations required for Batch Deconvolution %i\n', ite)])
LineBreakLarge()

if(SETTINGS.PLOT_GRAPHS)
    figure_name = sprintf('%s : Condition',mfilename);
    figure('name',figure_name)
    hold on
    plot(log10(condition),'-s')
    hold off
end




end

function Y = BuildY(arr_hx, vec_m)
% Build the coefficient matrix Y. This is the change of variable such
% that
% E(z)*h = Y(h)*f
%
% % Inputs
%
% arr_hx : (Array of Vectors) Vectors contain coefficients of the
% polynomials h_{i}(x)
%
% vec_m : (Vector) Degree of each polynomial f_{i}(x)
%

for i = 1:1:length(arr_hx)
    
    % Start with f1*h1
    % h_{1} is the first in the cell array h_{i}
    % f_{1} is the second in the cell array f_{i}
    % deg(f_{1}) = m(2)
    
    % Get polynomial h_{i}(x)
    hx = arr_hx{i};
    
    % Get degree of f_{i}
    deg_fw = vec_m(i+1);
    
    y{i} = real(BuildY1(hx, deg_fw));
end

% Build the Coefficient Matrix C

% Get number of columns for zero segment
nColumns = (vec_m(1)+1);


Y = blkdiag( y{1:length(y)});
nRows = size(Y,1);

Y = [zeros(nRows,nColumns) Y];




end

function Y1 = BuildY1(hx,m1)

% Construct a partition Y1 of the matrix Y.
Y1 = BuildT1(hx,m1);

end

function LHS_Matrix = BuildC(arr_fx,vMult)
% Build the LHS Coefficient matrix
%
% % Inputs
%
%
% arr_fx : (Array of Vectors) Array containing coefficients of polynomials
% f_{i}(x)
%
% vMult : (Vector) Contains multiplicity structure of the factors of f_{0}(x)



% For each distinct hx, build the partition of the Matrix C(f_{m_{i}+1},...,f_{m_{i+1}})
% C(f_{1},...f_{m1}), C(f_{m1+1},...,f_{m2}),
% where
% C(f_{1},...,f_{m1}) = [T(f_{1}) ; T(f_{2}) ; ... T(f_{m1})]

% Get number of distinct polynomials in h_{i}.
nDistinct_Polys_arr_hx = length(vMult);

arr_Cf = cell(nDistinct_Polys_arr_hx,1);

for i = 1 : 1 : nDistinct_Polys_arr_hx
    
    if i > 1
        old_mult = vMult(i-1);
    else % No previous multiplicity
        old_mult = 0;
    end
    
    new_mult = vMult(i);
    
    arr_Cf{i} = [];
    
    % for each polynomial f_{i} in the interval f_{m_{i-1}+1}...f_{m_{i}}
    arr_Tf = cell((new_mult+1) - (old_mult+1+1));
    
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
        
        % Stack beneath all other T_{f} which are multiplied by [_{i}(x)
        arr_Cf{i,1} = [arr_Cf{i,1} ; arr_Tf{j,1}];
    end
    
    
end

LHS_Matrix = blkdiag(arr_Cf{:});


end




function arr_hx = Get_hx(arr_px, vUniqueMult)

% Get number of entries in the array of polynomials p_{i}(x)
nEntries_arr_px = size(arr_px,1);

% initialise count
count = 1;


for i = 1:1:nEntries_arr_px
    
    if i == 1
        nReps = vUniqueMult(i);
    else
        nReps = (vUniqueMult(i) - vUniqueMult(i-1));
    end
    
    for j = 1:1:nReps
        arr_hx{count,1} = arr_px{i};
        count = count + 1;
    end
    
end
end


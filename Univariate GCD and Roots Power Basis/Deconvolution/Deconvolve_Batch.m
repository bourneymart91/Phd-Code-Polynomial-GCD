function arr_hx = Deconvolve_Batch(arr_fx)
% DECONVOLVE_BATCH Given the set of polynomials f_{0},...,f_{1}. compute
% the series of deconvolutions h_{1} = f_{1}/f_{0} h_{2} = f_{2}/f_{1},...
% Perform the deconvolutions by using the structure
% diag [ C(f_{1}) C(f_{2}) ... ] [h1 h2 ...]^{T} = [f_{0} f_{1} ...]^{T}
%
% % Inputs
%
% arr_fx : Array of polynomials f_{i}(x)
%
% % Outputs.
%
% arr_hx : Array of polynomials h_{i}(x)

global SETTINGS

% Get the number of polynomials in the set arr_fx
nPolys_arr_fx = size(arr_fx,1);
nPolys_arr_hx = size(arr_fx,1) - 1;

% Get the degree m_{i} of each of the polynomials f_{i} and store in a
% vector.
vDeg_arr_fx = zeros(nPolys_arr_fx, 1);
for i = 1:1:nPolys_arr_fx
    vDeg_arr_fx(i) = GetDegree(arr_fx{i});
end

% Get the degrees n{i} of polynomials h_{i} = f_{i}/f_{i+1}.
vDeg_arr_hx = zeros(nPolys_arr_fx-1, 1);
for i = 1:1:nPolys_arr_hx
    vDeg_arr_hx(i) = vDeg_arr_fx(i) - vDeg_arr_fx(i+1);
end

% Obtain theta such that the ratio of max element to min element is
% minimised

if( SETTINGS.PREPROC_DECONVOLUTIONS)
    theta = GetOptimalTheta(arr_fx,vDeg_arr_fx);
else
    theta = 1;
end

% Get polynomials f_{i}(\omega) by preprocessing
arr_fw = GetPolynomialArrayWithThetas(arr_fx,theta);

% Build RHS vector
RHS_vec = real(BuildRHSF(arr_fw));

% Build convolution matrix
DCQ = BuildC(arr_fw);

% Solve for initial values of h
hw_vec = SolveAx_b(DCQ,RHS_vec);

% Get array of polynomials given a vector of coefficients, and vector
% containing the degree of each polynomial.
arr_hw = GetPolyArrayFromVector(hw_vec, vDeg_arr_hx);

% Remove thetas to obtain array of polynomials h_{i}(x)
arr_hx = GetPolynomialArrayWithoutThetas(arr_hw, theta);


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
% arr_fx = array of polynomials f(x)
%
% % Outputs.
%
% C : Matrix DCQ containing coefficients of f_{i} as described above.

% Get the number of polynomials in the array f_{i}
nPolys_arr_fx = length(arr_fx);

% Initialise an array of cells to store the submatrices which will make up
% DCQ
T = cell(nPolys_arr_fx-1,1);

% For each polynomial f_{i}
for i = 1 : 1 : nPolys_arr_fx-1
    
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
    T{i} = BuildT1(fw, deg_hw);
end

% Build the Coefficient Matrix C of all matrices c
C = blkdiag(T{1:length(T)});

end



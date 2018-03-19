function arr_hx = Deconvolve_Batch(arr_fx)
% Perform batch deconvolution where each h_{i}(x) = f_{i-1}(x)/f_{i}(x) so
% each polynomial f_{i}(x) appears in the RHS matrix and the LHS vector,
% except f_{0} which only appears in the matrix and f_{n} which only
% appears in the vector.
%
% % Inputs
%
% arr_fx : (Array of Vectors) Array of polynomials f_{i}(x) to be used in a
% sequence of deconvolutios
%
% % Outputs
%
% arr_hx : (Array of Vectors) Array of polynomials h_{i}(x), the outputs of
% the deconvolutions h_{i} = f_{i}(x)/f_{i+1}(x)

global SETTINGS

% Get the number of polynomials f_{i}(x) in the set arr_fx
nPolynomials_arr_fx = size(arr_fx, 1);

% Intialise vector to store the degree of the set of polynomaisl f_{i}
vDegree_arr_fx = zeros(1,nPolynomials_arr_fx);

% For each polynomial f_{i}(x), get its degree and store in a vector
for i = 1 : 1 : nPolynomials_arr_fx
    
    vDegree_arr_fx(i) = GetDegree(arr_fx{i});
    
end

% Get the degrees n{i} of each polynomials h_{i} = f_{i}/f_{i+1}.
vDegree_arr_hx = (vDegree_arr_fx(1 : end - 1) - vDegree_arr_fx(2 : end))';


% Obtain theta such that the ratio of max element to min element is
% minimised
if(SETTINGS.PREPROC_DECONVOLUTIONS)
    
    theta = GetOptimalTheta(arr_fx);
    
    
else
    theta = 1;
end

% Get array of preprocessed polynomials f_{i}(\omega)
arr_fw = GetPolynomialArrayWithThetas(arr_fx, theta);

% Write Deconvolutions in form D^{-1}_{m+n} T_{m}(f(x))Q_{m} h = RHS_f

% Build the left hand matrix 
RHS_fw = BuildRHSF(arr_fw);

% Build the right hand side vector
DCQ = BuildDCQ(arr_fw);

% Get the solution vector h(w) in the system of equations
% DCQ * h = RHS_vec.
v_hw = SolveAx_b(DCQ, RHS_fw);


% Separate solution vector h, into component parts h_{1},h_{2},...h_{d},
% each of degree n_{i}
% initialise a cell array to store the coefficients of the individual
% polynomials h_{i}

% Split vec h in to an array of polynomials.
arr_hw = GetPolynomialArrayFromVector(v_hw, vDegree_arr_hx);

% Get polynomials h_{i}(x)
arr_hx = GetPolynomialArrayWithoutThetas(arr_hw, theta);



end




function DCQ = BuildDCQ(arr_fx)
% set fw is the cell array of poly coefficiencts fw_i
%
% Inputs.
%
% arr_fx : (Array of Vectors) Each cell contains coefficients of polynomial
% f_{i}(x)
%
% % Outputs
%
% DCQ : (Matrix) Convolution matrix in the deconvolution problem

% Get the number of polynomials in the array
nPolynomials_arr_fx = length(arr_fx);

% Initialise the cell array to store the matrices D^{-1}_{}T_{}(f)Q_{}
arr_DT1Q1 = cell(nPolynomials_arr_fx - 1,1);

% For each of the polynomials f_{i}(x), excluding the final polynomial
for i = 1 : 1 : (nPolynomials_arr_fx - 1)
    
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
    
    % Build the matrix D^{-1}
    D = BuildD_2Polys(deg_fw, deg_hw);
    
    % Build the Matrix T_{}(f)
    T1 = BuildT1(fw, deg_hw);
    
    % Build the matrix Q_{}
    Q1 = BuildQ1(deg_hw);
    
    % Add matrix DTQ to array of matrices.
    arr_DT1Q1{i}  = D*T1*Q1;
    
end

% Build the coefficient matrix DTQ
DCQ = blkdiag(arr_DT1Q1{1:length(arr_DT1Q1)});

end


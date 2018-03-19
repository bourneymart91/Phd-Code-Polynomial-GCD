function arr_hxy = Deconvolve_Bivariate_Batch_Both(arr_fxy, vDeg_t_arr_fxy, vDeg_x_arr_fxy, vDeg_y_arr_fxy)
% Perform the batch deconvolution with the constraints on the total degree
% and relative degree
%
% Inputs.
%
% arr_fxy : (Array) Array of polynomials f(x,y)
%
% vDegt_arr_fxy : (Vector) Vector containing the total degree of the
% polynomials f_{i}(x,y)
%
% vDegx_arr_fxy : (Vector) Vector containing the degree of polynomials
% f_{i}(x,y) with respect to x
%
% vDegy_arr_fxy : (Vector) Vector containing the degree of polynomials
% f_{i}(x,y) with respect to y
%
% % Outputs
%
% arr_hxy : (Array) Array of polynomials h(x,y)

% Form the left hand side matrix

% Get the number of polynomials in the array arr_f_{i}(x,y)
nPolys_arr_fxy = size(arr_fxy,1);
nPolys_arr_hxy = nPolys_arr_fxy - 1;

% Get the total degree of polynomials in the array h_{i}(x,y)
vDeg_t_arr_hxy = abs(diff(vDeg_t_arr_fxy));

% Get the degree of polynomials h_{i}(x,y) with respect to x and y
vDeg_x_arr_hxy = abs(diff(vDeg_x_arr_fxy));
vDeg_y_arr_hxy = abs(diff(vDeg_y_arr_fxy));

% Initialise an array to store convolution matrices 
arr_T1 = cell(nPolys_arr_hxy, 1);

% For each of the polynomials in the array f_{i}(x,y) excluding the first,
% build the convolution matrix
for i = 2 : 1 : nPolys_arr_fxy
    
    fxy = arr_fxy{i};
    m_current = vDeg_t_arr_fxy(i);
    
    n  = vDeg_t_arr_hxy(i-1);
    n1 = vDeg_x_arr_hxy(i-1);
    n2 = vDeg_y_arr_hxy(i-1);
    
    % Build the matrix T(f(x,y))
    T1 = BuildT1_Both_Bivariate(fxy, m_current, n, n1, n2);
    
    % Add T1 to array of matrices
    arr_T1{i-1} = T1;
    
end

% Build a block diagonal of the matrices
C_fww = blkdiag(arr_T1{:});

vRHS = BuildRHS_Vector(arr_fxy, vDeg_t_arr_fxy, vDeg_x_arr_fxy, vDeg_y_arr_fxy);


% %
x_ls = pinv(C_fww) * vRHS;
%x_ls = SolveAx_b(C_fww,vRHS);


% Initialise array to store polynomials h_{i}(x,y)
arr_hxy = cell(nPolys_arr_hxy, 1);

for i = 1 : 1 : nPolys_arr_hxy
    
    
    % Ge the total degree of h_{i}(x,y)
    n = vDeg_t_arr_hxy(i);
    
    % Get degree of h_{i}(x,y) with respect to x and y
    n1 = vDeg_x_arr_hxy(i);
    n2 = vDeg_y_arr_hxy(i);
    
    % Get number of nonzeros 
    nNonZeros_hxy = GetNumNonZeros(n1,n2,n);
    
    % Get the coefficients of h_{i}(x,y)
    temp_vec = x_ls(1:nNonZeros_hxy);
    
    % Remove the subset of coefficients from the vector x_ls
    x_ls(1:nNonZeros_hxy) = [];
    
    % Append zeros so that temp_vec fills a n1+1,n2+1 matrix
    nCoefficients = (n1+1) * (n2+1);
    nZeros = nCoefficients - nNonZeros_hxy;
    vCoefficients = [temp_vec ; zeros(nZeros,1)];
    
    arr_hxy{i} = GetAsMatrix(vCoefficients, n1, n2);
    
    
    
    
end
end


function vRHS = BuildRHS_Vector(arr_fxy, vDeg_t_arr_fxy, vDeg_x_arr_fxy, vDeg_y_arr_fxy)
% Form the Right hand side vector
% 
% % Inputs
%
% arr_fxy : (Array of Matrices) Each cell in the array contains
% coefficients of a polynomial f_{i}(x,y)
%
% vDeg_t_arr_fxy : (Vector) Contains total degree of each polynomial
% f_{i}(x,y)
%
% vDeg_x_arr_fxy : (Vector) Contains degree of each polynomial f_{i}(x,y)
% with respect to x
%
% vDeg_y_arr_fxy : (Vector) Contains degree of each polynomial f_{i}(x,y)
% with respect to y
%
% % Outputs
%
% vRHS : (Vector) Right hand side vector containing coefficients of
% f_{0},..., f_{n-1}

% Get the number of polynomials in the array
nPolys_arr_fxy = length(arr_fxy);

% Initialise an array to store coefficients of the polynomials f_{i}(x,y)
% as vectors.
arr_rhs = cell(nPolys_arr_fxy - 1,1);

for i = 1:1:nPolys_arr_fxy - 1
    
    % temporarily label the ith entry of the array as fxy
    v_fxy = GetAsVector(arr_fxy{i});
    
    % Get total and relative degree structure of f(x,y)
    m = vDeg_t_arr_fxy(i);
    m1 = vDeg_x_arr_fxy(i);
    m2 = vDeg_y_arr_fxy(i);
    
    % Remove zeros from the end of the vector
    nNonZeros_fxy = GetNumNonZeros(m1, m2, m);
    
    % Remove the zero entries
    v_fxy = v_fxy(1:nNonZeros_fxy,1);
    
    % Add the vector to the array
    arr_rhs{i} = v_fxy;
end

% Get RHS Vector by converting cell array to matrix
vRHS = cell2mat(arr_rhs);

end
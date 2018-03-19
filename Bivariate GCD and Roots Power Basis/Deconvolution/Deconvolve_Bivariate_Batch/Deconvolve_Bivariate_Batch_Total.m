function arr_hxy = Deconvolve_Bivariate_Batch_Total(arr_fxy,vDegt_arr_fxy)
% Perform the batch deconvolution C(f1,f2,..,fn-1) h = [f0;f1;,...,;fn]
%
% Inputs.
%
% arr_fxy : (Array of Matrices) Array of polynomials f(x,y)
%
% vDeg_fxy : (Vector) Vector containing the total degree of the polynomials
% f{i}

% %
% %
% Form the left hand side matrix

% Get total degree of the polynomials h_{i}(x,y)
vDegt_arr_hxy = abs(diff(vDegt_arr_fxy));


% % Build the convolution matrix
C = BuildConvolutionMatrix(arr_fxy, vDegt_arr_fxy);

% % Build the right hand side vector
vRHS = BuildRHS_Vector(arr_fxy, vDegt_arr_fxy);

% Solve Ax = b
vec_hxy = SolveAx_b(C,vRHS);

% Initialise cell array 
arr_hxy = GetPolynomialArrayFromVector(vec_hxy, vDegt_arr_hxy);

end


function C = BuildConvolutionMatrix(arr_fxy, vDegt_arr_fxy)
% 
%
% % Inputs
%
% arr_fxy : (Array of Matrices) 
%
% vDegt_arr_fxy : (Vector) Total degree of each polynomial f_{i}(x,y)

% Get number of polynomials in the array
nPolys_arr_fxy = length(arr_fxy);

% Get the degree of polynomials h_{i}(x,y)
vDegt_arr_hxy = abs(diff(vDegt_arr_fxy));


% Initialise cell array to store matrices T_{}(f_{i}(x,y))
arr_T1 = cell(nPolys_arr_fxy - 1,1);

% For each of the polynomials excluding the first.
for i = 2:1:nPolys_arr_fxy 
    
     
    % Build the matrix T(f(x,y))
    T1 = BuildT1_Total_Bivariate(arr_fxy{i,1},vDegt_arr_fxy(i,1),vDegt_arr_hxy(i-1,1));
    
    % % Strip the columns corresponding to zeros in f{i-1}
    arr_T1{i-1} = T1;
end

C = blkdiag(arr_T1{:});

end

function vRHS = BuildRHS_Vector(arr_fxy, vDegt_arr_fxy)
%
% % Inputs
%
% arr_fxy : (Array of Matrices) Each cell contains a matrix of coefficients
% of a polynomial f_{i}(x,y) for i = 0,...,n
% 
% % Outputs
%
% vRHS : (Vector) Contains coefficients of polynomials f_{0},...,f_{n-1}


% Get number of polynomials in the array
nPolys_arr_fxy = length(arr_fxy);

% Initialise an array to store vectors of coefficients
arr_rhs = cell(nPolys_arr_fxy -1,1);

% Form the Right hand side vector
for i = 1:1:nPolys_arr_fxy - 1
    
    % Temporarily label the ith entry of the array as fxy
    v_fxy = GetAsVector_Version1(arr_fxy{i});
    
    % Remove the corresponding zeros
    
    % Get total degree of polynomial
    m = vDegt_arr_fxy(i);
    
    % Get number of nonzero coefficients
    nNonZeros_fxy = nchoosek(m+2,2);
    
    % Remove zeros from vector
    v_fxy = v_fxy(1:nNonZeros_fxy);
    
    arr_rhs{i} = v_fxy;
    
end

vRHS = cell2mat(arr_rhs);

end

function arr_hxy = GetPolynomialArrayFromVector(vec_hxy, vDegt_arr_hxy)
%
%
% % Inputs
%
% vec_hxy : (Vector) Coefficients of all polynomials h_{i}(x,y)
%
% vDegt_arr_hxy : (Vector) Total degree of each of the polynomials
% h_{i}(x,y)
%
% % Outputs
%
% arr_hxy : (Array of Vectors) Each cell in the array contains matrix of
% coefficients of a polynomial h_{i}(x,y)

% Get number of polynomials in the array
nPolys_arr_hxy = length(vDegt_arr_hxy);

% Initialise the cell array 
arr_hxy = cell(nPolys_arr_hxy,1);

for i = 1:nPolys_arr_hxy
    
    
    n = vDegt_arr_hxy(i);
    nNonZeros = nchoosek(n+2,2);
    
    % Get vector of zeros 
    temp_vec = zeros((n+1)*(n+1),1);
    temp_vec(1:nNonZeros) = vec_hxy(1:nNonZeros);
    
    mat = GetAsMatrix_Version1(temp_vec, n, n);
    
    arr_hxy{i} = mat;
    vec_hxy(1:nNonZeros) = [];
    
end

end

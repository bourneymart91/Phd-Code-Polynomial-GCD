function arr_hxy = Deconvolve_Bivariate_Batch_Respective(arr_fxy)
% Perform the batch deconvolution
%
% Inputs.
%
% arr_fxy : (Array of Matrices) Array of polynomials f(x,y)
%
%
% % Outputs
%
% arr_hxy : (Array of Matrices) Each cell contains a matrix of coefficients
% of polynomial f(x,y)



% Form the left hand side matrix

% Get number of polynomials in the array arr_fxy
nPolys_arr_fxy = size(arr_fxy,1);

% For each of the polynomials excluding the first.
vNCoefficients_hxy = zeros(nPolys_arr_fxy-1,1);

% Initialise vectors to store degree of polynomials f_{i}(x,y) with respect
% to x and y
vDeg_x_fxy = zeros(nPolys_arr_fxy, 1);
vDeg_y_fxy = zeros(nPolys_arr_fxy, 1);

for i = 1:1:nPolys_arr_fxy

    % Get degree of each polynomial in the array
    [vDeg_x_fxy(i), vDeg_y_fxy(i)] = GetDegree_Bivariate(arr_fxy{i});
    
end

nPolys_arr_hxy = nPolys_arr_fxy - 1;

% Get degree structure of polynomials in array h_{i}(x,y)
vDeg_x_hxy = abs(diff(vDeg_x_fxy));
vDeg_y_hxy = abs(diff(vDeg_y_fxy));

arr_T1 = cell(nPolys_arr_hxy,1);

for i = 2:1:nPolys_arr_fxy
    
   
    % Temporarily call the ith entry f(x,y)
    fxy = arr_fxy{i};
    

    % Get the degree of h{i-1}(x,y)
    n1 = vDeg_x_hxy(i-1);
    n2 = vDeg_y_hxy(i-1);
    
    % Get number of coefficients in h_{i}(x,y)
    vNCoefficients_hxy(i-1) = (n1 + 1) * (n2+1);
    
    % Build the matrix T(f(x,y))
    T1 = BuildT1_Relative_Bivariate(fxy, n1, n2);
     
    arr_T1{i - 1} = T1;
end

% Build the coefficient matrix
C = blkdiag(arr_T1{:});

% %
% %
% Form the Right hand side vector
vRHS = GetRHS_Vector(arr_fxy);

% %
% %
% %
x_ls = SolveAx_b(C,vRHS);

% Split solution vector into polynomials h{i}.

arr_hxy = cell(nPolys_arr_hxy,1);

for i = 1:nPolys_arr_hxy
    
    temp_vec = x_ls(1:vNCoefficients_hxy(i));
    
    arr_hxy{i,1} = GetAsMatrix(temp_vec, vDeg_x_hxy(i), vDeg_y_hxy(i));
    
    x_ls(1:vNCoefficients_hxy(i)) = [];    
end

end


function vRHS = GetRHS_Vector(arr_fxy)
%
% % Inputs
%
% arr_fxy : (Array of Matrices) Each cell contains matrix of coefficients
% of polynomail f_{i}(x,y)
%
% % Outputs
%
% vRHS : (Vector) Vector containing coefficients of all polynomials
% f_{0}(x,y),...,f_{n-1}(x,y)


% Get number of polynomials in the array 
nPolys_arr_fxy = length(arr_fxy);


% Initialise cell array to store coefficients of f_{0},...,f_{d-1}
arr_rhs = cell(nPolys_arr_fxy - 1,1);

for i = 1:1:nPolys_arr_fxy - 1
    
    % Add polynomial to array
    arr_rhs{i} = GetAsVector(arr_fxy{i});
    
end

vRHS = cell2mat(arr_rhs);


end
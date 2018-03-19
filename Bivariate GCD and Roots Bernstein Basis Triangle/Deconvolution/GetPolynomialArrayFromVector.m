
function arr_hxy = GetPolynomialArrayFromVector(v_hxy, vDeg_arr_hxy)
%
% % Inputs
%
% v_hxy :
%
% vDeg_arr_hxy

% Get number of polynomials in array v_{x,y}
nPolys_arr_hxy = size(vDeg_arr_hxy,1);

% Initialise a cell array to store h_{i}(x,y).
arr_hxy = cell(nPolys_arr_hxy,1);


for i = 1:1:nPolys_arr_hxy
    
    % the polynomial h_{i} has degree m_{i-1} - m_{i}
    n = vDeg_arr_hxy(i);
    
    % Number of coefficients in h_{i}(x,y)
    nCoefficients_hxy = nchoosek(n+2,2);
    
    % Vector of coefficients of h_{i}(x,y)
    vec = v_hxy(1:nCoefficients_hxy);
    
    % append zeros to vector so we can form a matrix
    try
        nZeros_hxy = nchoosek(n+1,2);
    catch
        nZeros_hxy = 0;
    end
    
    
    vec = ...
        [
        vec;
        zeros(nZeros_hxy,1)
        ];
    
    % Form matrix of coefficients of h_{i}(x,y)
    arr_hxy{i} = GetAsMatrix(vec,n,n);
    
    % Remove coefficients from solution vector x_ls
    v_hxy(1:nCoefficients_hxy) = [];
end
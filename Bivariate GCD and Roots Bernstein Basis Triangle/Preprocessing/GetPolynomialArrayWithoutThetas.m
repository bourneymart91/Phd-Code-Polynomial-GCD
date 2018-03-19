function arr_fxy = GetPolynomialArrayWithoutThetas(arr_fww, th1, th2)
%
% % Inputs
% 
% arr_hww : (Array of Matrices) Each cell in the array contains matrix of
% coefficients of the polynomial f_{i,j}(\omega_{1},\omega_{2})
%
% th1 : (Float)
%
% th2 : (Float)
%
% % Outputs
%
% arr_hxy : (Array of Matrices) Each cell in the array contains matrix of
% coefficients of the polynomial f_{i,j}(x,y)


% Get number of polynomials in the array
nPolys_arr_fxy = length(arr_fww);

% Initialise cell array
arr_fxy = cell(nPolys_arr_fxy, 1);

for i = 1:1:nPolys_arr_fxy
    
    [m,~] = GetDegree_Bivariate(arr_fww{i});
    
    arr_fxy{i} = GetWithoutThetas(arr_fww{i}, m, th1, th2);
    
end


end
function arr_hx = GetPolynomialArrayWithoutThetas(arr_hw, theta)
% Get array of polynomials h_{i}(x) from h_{i}(\omega)
%
% % Inputs
%
%
% arr_hw : (Array of Vectors) Array of polynomials h_{i}(\omega), where
%   the i-th cell contains coefficients of h_{i}(\omega)
%
% theta : (Float) Value of \theta
%
% % Outputs
%
%
% arr_hx : (Array of Vectors) Array of polynomials h_{i}(x) where the ith 
%   cell contains coefficients of h_{i}(x)



% Get the number of polynomials in the array
nPolynomials_arr_hx = length(arr_hw);


% Initialise a cell array to store the polynomials h_{i}(x)
arr_hx = cell(nPolynomials_arr_hx, 1);

% For each polynomial in the array h_{i}(\omega)
for i = 1:1:nPolynomials_arr_hx
    
    % Remove thetas and store in array of polynomials h_{i}(x)
    arr_hx{i} = GetWithoutThetas(arr_hw{i}, theta);
    
end

end
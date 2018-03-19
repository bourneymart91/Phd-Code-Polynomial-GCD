function arr_hx = GetPolynomialArrayWithoutThetas(arr_hw, theta)
% Given an array of polynomials h_{i}(\omega), get h_{i}(x)
%
% % Inputs
%
% arr_hw : (Array of Vectors) Array of vectors containing coefficients of
% polynomials h_{i}(\omega)
%
% theta : (Float) Optimal value of \theta
%
% % Outputs
%
% arr_hx : (Array of Vectors) Array of vectors containing coefficietns of
% the polynomials h_{i}(x)

% Get number of polynomials in the array
nPolys_arr_hx = length(arr_hw);

% Remove thetas from h_{i}(w) to get h_{i}(x)
arr_hx = cell(nPolys_arr_hx,1);

for i = 1:1:nPolys_arr_hx

    arr_hx{i} = GetWithoutThetas(arr_hw{i}, theta);

end

end


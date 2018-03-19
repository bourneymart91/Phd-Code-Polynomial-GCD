function arr_fw = GetPolynomialArrayWithThetas(arr_fx, theta)
% Given an array of polynomials f_{i}(x) for i = 0,...,n. Get each
% polynomial in preprocessed form f_{i}(\omega).
%
% % Inputs
%
% arr_fx : (Array of Vectors) Coefficients of polynomials f_{i}(x)
%
% theta : (Float) \theta
%
% % Outputs
%
% arr_fw : (Array of Vectors) Coefficients of polynomials f_{i}(\omega)

nPolys_arr_fx = length(arr_fx);

% Initialise a cell-array for f(w)
arr_fw = cell(nPolys_arr_fx, 1);

% for each f_{i} get fw_{i}
for i = 1 : 1 : nPolys_arr_fx
    arr_fw{i} = GetWithThetas(arr_fx{i}, theta);
end


end
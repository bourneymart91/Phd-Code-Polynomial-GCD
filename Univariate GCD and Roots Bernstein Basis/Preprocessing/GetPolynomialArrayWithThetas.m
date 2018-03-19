function arr_fw = GetPolynomialArrayWithThetas(arr_fx, theta)
% Get an array of polynomials f_{i}(x) in the form f_{i}(\omega) 
%
% % Inputs
%
% arr_fx : (Array of Vectors) Each cell of the array contains coefficients 
%   of the polynomial f_{i}(x).
%
% theta : (Float)
%
% % Outputs
%
% arr_fw : (Array of Vectors) Each cell of the array contains coefficients
%   of the preprocessed polynomial f_{i}(\omega).

% Get the number of polynomials in the array
nPolys_arr_fx = length(arr_fx);

% Initialise a cell-array for f(\omega) the prepfocessed form of f(x)
arr_fw = cell(nPolys_arr_fx, 1);

% For each polynomial f_{i}(x), preprocess to obtain f_{i}(\omega)
for i = 1 : 1 : nPolys_arr_fx
    
    arr_fw{i,1} = GetWithThetas(arr_fx{i}, theta);
    
end


end
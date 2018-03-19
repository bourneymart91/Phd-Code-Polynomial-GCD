function vDegree = GetDegree_Array(arr_fx)
% Get the degree of an array of polynomials f_{i}(x)
% 
% % Inputs
%
%
% arr_fx : (Array of Vectors) An array of vectors where each cell in the
% array contains the vector of the coefficients of f_{i}(x)
%
%
% % Outputs
%
% vDegree : (Vector) A vector whose i-th entry contains the degree of the
% ith polynomial in the array arr_fx


% Get number of polynomials in the array
nPolys = length(arr_fx);

% Initialise a vector to store the degrees of each polynomial
vDegree = zeros(nPolys, 1);

% For each polynomial in the array get its degree
for i = 1 : 1: nPolys
    
    fx = arr_fx{i};
    
    vDegree(i) = GetDegree(fx);
    
end



end
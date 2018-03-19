function arr_fx = GetPolyArrayFromVector(vec_fx, vDegree_fx)
% Get array of polynomials f_{i}(x), given a vector of the coefficients
% containing all polynomials, and the degree of each polynomial. THe vector
% must contain coefficients in the order [f_{0} f_{1} ...]
%
% % Inputs
%
% vec_fx : (Vector) Contains coefficients of polynomials
% f_{0}(x),...,f_{n-1}(x)
%
% vDegree_fx : (Vector) Contains degree of each polynomial f(x)
%
%
% % Outputs
%
% arr_fx : (Array of Vectors) Array where each cell contains coefficients
% of each polynomial f_{i}(x)

% Given the vector of perturbations of f(x) given by v_zx

% Get number of polynomials in arr_fx
nPolys_arr_fx = length(vDegree_fx);

% Initialise an array
arr_fx = cell(nPolys_arr_fx, 1);

for i = 1:1:nPolys_arr_fx
    
    % Get the m+1 coefficients from the vector
    arr_fx{i} = vec_fx(1:vDegree_fx(i)+1);
    
    % Remove the m+1 coefficients
    vec_fx(1:vDegree_fx(i)+1) = [];
    
end


end
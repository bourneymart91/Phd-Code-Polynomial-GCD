function arr_fx = GetPolynomialArrayFromVector(vec_fx, vDegree_fx)
% Given the vector of perturbations of f(x) given by v_zx
%
% % Inputs
%
% vec_fx : (Vector) Contains coefficients of all polynomials f_{0},...,f_{n}
%
% vDegree_fx : (Vector) Contains degree of each polynomial f_{0},...,f_{n}
%
% % Outputs
%
% arr_fx : (Array of Vectors) Each cell of the array contains coefficients
% of the polynomial f_{i}(x).


% Get number of polynomials in arr_fx
nPolys_fx = length(vDegree_fx);

% Initialise an array
arr_fx = cell(nPolys_fx, 1);

for i = 1:1:nPolys_fx
    
    % Get the m+1 coefficients from the vector
    arr_fx{i} = vec_fx(1 : vDegree_fx(i) + 1);
    
    % Remove the m+1 coefficients
    vec_fx(1:vDegree_fx(i)+1) = [];
    
end


end
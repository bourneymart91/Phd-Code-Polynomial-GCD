
function f = BuildRHSF(arr_fw)
% Build the vector f such that it contains the elements of
% Rhs f = [f_{0},...,f_{n-1}]. This vector forms part of the deconvolution
% problem C(f) h = f
%
% % Inputs
%
% fw = (Array of Vectors) array of vectors f_{0},...,f_{n}
%
% % Outputs
%
% f : (vector) Coefficients of the polynomials f_{0},..., f_{n-1}

% Initialise empty vector.
f = [];

% Get number of polynomials in the array f_{i}(x)
nPolys_arr_fw = length(arr_fw);

% Add all but the last polynomial to a vector
for i = 1 : 1 : (nPolys_arr_fw - 1)
    
    f = [ f; arr_fw{i}];
    
end

end
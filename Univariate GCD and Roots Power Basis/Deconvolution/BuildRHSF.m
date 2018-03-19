function f = BuildRHSF(arr_fx)
% Build the vector f such that it contains the elements of
% Rhs f = [f_{0},...,f_{n-1}]
%
%
% fw = array of vectors f_{0},...,f_{n}
%

% Initialise empty vector.
f = [];

% Get number of polynomials in the array
nPolys_arr_fx = length(arr_fx);


% for each vector f f_{0},...,f_{n-1} in fw_array, add to right hand
% side vector

for i= 1 : 1 : nPolys_arr_fx - 1
    f = [f;arr_fx{i}];
end

end
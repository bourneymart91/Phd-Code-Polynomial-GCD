function arr_hx = Deconvolve_Separate(arr_fx)
% Given a set of polynomials f_{i}, Perform a series of deconvolutions
% h_{i}(x) = f_{i}(x)/f_{i+1}(x)
%
% Each deconvolution is performed independently, and the solutions are
% output as an array of vectors h_{i}(x).
%
% Inputs
%
%
% arr_fx : (Array of Vectors) Array of Vectors containing coefficients of
% polynomials.
%
% % Outputs
%
% arr_hx : (Array of Vectors) Array of vectors containing coefficients of
% polynomials h_{i}(x) where each h_{i}(x) = f_{i}(x)/f_{i+1}(x)


% Get number of polynomials in array of f_{i}(x)
nPolys_arr_fx = size(arr_fx, 1);

% Initialise the array of h_{i}(x)
arr_hx = cell(nPolys_arr_fx - 1, 1);

% for each item in set g starting at
for i = 1 : 1 : nPolys_arr_fx - 1
    
    % Get the two polynomials which are to be deconvolved
    % f{i},f{i+1}
    arr_hx{i,1} = Deconvolve(arr_fx{i}, arr_fx{i+1});
    
end


end

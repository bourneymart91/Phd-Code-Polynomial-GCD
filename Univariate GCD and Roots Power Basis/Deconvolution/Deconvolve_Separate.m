function arr_hx = Deconvolve_Separate(arr_fx)
% Peform the sequence of deconvolutions independently
%
% % Inputs
%
% arr_fx : (Array of Vectors) Vectors containing coefficients of
% polynomials f_{i}(x)
%
%
% % Outputs
%
% arr_hx : (Array of Vectors) Vectors containing coefficients of
% polynomials h_{i}(x)


% Get number of polynomials f_{i}(x)
nPolys_arr_fx = length(arr_fx);

% Initialise an array to store polynomials h_{i}(x) = h{i-1} ./ h{i}
arr_hx = cell(nPolys_arr_fx-1,1);

% For each pair of polynomials perform deconvolution.
for i = 1:1:nPolys_arr_fx-1
    
    arr_hx{i} = Deconvolve(arr_fx{i},arr_fx{i+1}) ;
    
end

end
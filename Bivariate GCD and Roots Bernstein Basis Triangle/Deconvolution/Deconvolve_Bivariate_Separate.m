function arr_hxy = Deconvolve_Bivariate_Separate(arr_fxy)
%
% % Inputs
%
% arr_fxy : (Array of Matrices) 
%
%
% % Outputs
%
% arr_hxy
%
%
%
% Peform the sequence of deconvolutions independently

% Get number of polynomials f_{i}(x)
nPolynomials_fxy = length(arr_fxy);

% Initialise an array to store polynomials h_{i}(x) = h{i-1} ./ h{i}
arr_hxy = cell(nPolynomials_fxy - 1, 1);

% For each pair of polynomials perform deconvolution.
for i = 1 : 1 : nPolynomials_fxy - 1
    
    arr_hxy{i,1} = Deconvolve_Bivariate(arr_fxy{i, 1}, arr_fxy{i + 1, 1}) ;
    
end

end
function arr_fww = GetPolynomialArrayWithThetas(arr_fxy, th1, th2)
% Get f(\omega_{1},\omega_{2}) from f(x,y)
%
% % Inputs
%
% arr_fxy : (Array of Matrices)
%
% th1 : (Float)
%
% th2 : (Float)
%
% % Outputs
%
% arr_fww : (Array of Matrices) 



% Get the number of polynomials in the array 
nPolys_arr_fxy = length(arr_fxy);


% Initialise a cell array to store preprocessed polynomials
arr_fww = cell(nPolys_arr_fxy,1);

for i = 1:1:nPolys_arr_fxy
    
    [m,~] = GetDegree_Bivariate(arr_fxy{i});
    
    arr_fww{i} = GetWithThetas(arr_fxy{i}, m, th1, th2);
    
end


end

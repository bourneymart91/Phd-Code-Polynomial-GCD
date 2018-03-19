
function vRHS = BuildRHS_vec(arr_fxy)
%
% % Inputs

% Get number of polynomials in array f_{i}(x,y)
nPolys_arr_fxy = size(arr_fxy,1);

% Get the degree of the polynomials f_{i}(x,y)
vDeg_arr_fxy = zeros(nPolys_arr_fxy,1);

for i = 1:1:nPolys_arr_fxy
    
    vDeg_arr_fxy(i) = GetDegree_Bivariate(arr_fxy{i});
    
end


% Initialise an array
arr_rhs = cell(nPolys_arr_fxy,1);

% For each polynomial f_{1},...,f_{d} (Note exclude f_{0})
for i = 1:1:nPolys_arr_fxy -1
    
    fxy = GetAsVector(arr_fxy{i});
    
    % Get degree of f(x,y)
    m = vDeg_arr_fxy(i);
    
    % Get number of coefficients in f(x,y)
    nCoefficients_fxy = nchoosek(m+2,2);
    
    % Strip zeros from f(x,y)
    v_fxy = fxy(1:nCoefficients_fxy);
    
    arr_rhs{i} = v_fxy;
    
end

vRHS = cell2mat(arr_rhs);

end
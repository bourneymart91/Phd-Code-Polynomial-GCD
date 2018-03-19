function [arr_hx] = GetArray_hx(fx_factor_multiplicity_matrix)
%
% % Inputs
%
% fx_factor_multiplicity_matrix : (Matrix) containing symbolic factors and
% their multiplicities
%


% Get multiplicity structure of f_{0}(x)
vMultiplicity_f0 = double(fx_factor_multiplicity_matrix(:,2));

% Get the set of factors of f_{0}(x)
vFactors = (fx_factor_multiplicity_matrix(:,1));

% Get multiplicity structure of f_{i}(x)
arr_vMultiplicity_fx = GetMultiplicityArr_fx(vMultiplicity_f0);

nPolysArr_fx = max(vMultiplicity_f0) + 1;

% Get number of polynomials h_{i}(x)
nPolys_arr_hx = nPolysArr_fx - 1;

% Initialise array to store h_{i}
arr_hx = cell(nPolys_arr_hx, 1);

for i = 1 : 1 : nPolys_arr_hx
   
    vMultiplicity_hi = arr_vMultiplicity_fx{i} - arr_vMultiplicity_fx{i+1};
    
    
    arr_hx{i} = BuildPolyFromRootsSymbolic([vFactors vMultiplicity_hi]);
    
    
end



end

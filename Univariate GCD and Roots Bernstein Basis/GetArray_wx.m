function [arr_wx] = GetArray_wx(fx_factor_multiplicity_matrix)
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

% Get number of polynomials w_{i}(x)
nPolys_arr_wx = nPolys_arr_hx;


% Initialise array to store h_{i}
arr_wx = cell(nPolys_arr_hx, 1);



% Get maximum multiplicity
maxMultiplicity = max(vMultiplicity_f0);

% Initialise a cell array, each cell contains a vector of symbolic factors.
arr_vFactors = cell(maxMultiplicity,1);

nFactors = size(fx_factor_multiplicity_matrix, 1);

for i = 1 : 1 : nFactors
   
    % Get factor
    myFactor = fx_factor_multiplicity_matrix(i,1);
    
    % Get multiplicity of factor
    myMultiplicity = fx_factor_multiplicity_matrix(i,2);
    
    % Get vector of factors of multiplicity 'myMultiplicity' 
    vFactors = arr_vFactors{myMultiplicity};
    vFactors = [vFactors ; myFactor];
    arr_vFactors{myMultiplicity} = vFactors;
    
end



for i = 1 : 1 : maxMultiplicity
   
    vFactors = arr_vFactors{i};
    vMultiplicity = ones(length(vFactors),1);
    
    if length(vFactors) >= 1
        arr_wx{i} = BuildPolyFromRootsSymbolic([vFactors vMultiplicity]);
    else
        arr_wx{i} = 1;
    end
end






end

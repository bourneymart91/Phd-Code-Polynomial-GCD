function [arr_fx] = GetArray_fx(fx_factor_multiplicity_matrix)
%
%
% fx_factor_multiplicity_matrix : Matrix containing the symbolic factors of
% f(x) and the multiplicity of the factors.
%
% % Outputs
%
% arr_fx : (Array of Vectors) Each vector contains the coefficients of a
% polynomial f_{i}(x) in the sequence generated in Gauss method of
% factorisation.
%
%

% Get symbolic factors
vFactors = (fx_factor_multiplicity_matrix(:,1));

% Get multiplicity of factors
vMultiplicity_f0 = double(fx_factor_multiplicity_matrix(:,2));

% Get max multiplicty
maxMultiplicity = max(vMultiplicity_f0);

% Set number of polynomials in array
nPolys_arr_fx = maxMultiplicity + 1;

% Initialise array
arrMultiplicities = cell(nPolys_arr_fx, 1);
arr_fx = cell(nPolys_arr_fx, 1);

for i = 1 : 1 : nPolys_arr_fx
    
    % Get multiplicity of factors in f_{i}(x)
    vMultiplicity_fi = vMultiplicity_f0 - (i-1);
    vMultiplicity_fi(vMultiplicity_fi < 0) = 0;
    
    arrMultiplicities{i} = vMultiplicity_fi;
    
    % Get coefficients of f_{i}(x)
    arr_fx{i} = BuildPolyFromRootsSymbolic([vFactors vMultiplicity_fi]);
    
end


end
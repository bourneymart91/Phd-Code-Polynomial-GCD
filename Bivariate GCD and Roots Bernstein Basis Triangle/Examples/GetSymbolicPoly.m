
function [symbolic_fxy] = GetSymbolicPoly(root_mult_arr)
% Given the set of symbolic factors of f(x,y) return the symbolic
% polynomial.


nUniqueFactors = size(root_mult_arr,1);

% Intialise the symbolic polynomial
symbolic_fxy = 1;

for i = 1:1:nUniqueFactors
    
    factor = root_mult_arr(i,1);
    mult = root_mult_arr(i,2);
    
    symbolic_fxy = symbolic_fxy * (factor.^mult);
end


end
function m = GetDegree_Univariate(fx)

% Ensure that fx is a column vector

[nRows,nCols] = size(fx);

if nCols ~= 1 
    error('Coefficients of polynomial are not a column vector');
end

m = nRows - 1;

end
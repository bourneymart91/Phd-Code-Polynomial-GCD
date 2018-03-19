function [lambda] = GetGeometricMean_Total(fxy,m)
%
% Inputs
%
% fxy : Matrix of coefficients of f(x,y
%

% Get f(x,y) as a vector
f = GetAsVector(fxy);

% Remove the zeros from the vector of coefficients of f(x,y)
nCoefficients_fxy = nchoosek(m+2,2);

f = f(1:nCoefficients_fxy);

% Get the number of zero coefficients in f(x,y)
nZeros = sum(f(:)==0);

% Get geometric mean of the non-zero entries
lambda = prod(fxy(fxy~=0)) .^(1./nchoosek(m+2,2)-nZeros);

%lambda = prod(fxy(fxy~=0).^ (1./nchoosek(m+2,2))) ;




end

function lambda = GetGeometricMean_Total_ByLogs(fxy,m)

% Remove zero values from f(x,y) 
fxy = abs(fxy(fxy~=0));

% Get f(x,y) in logs
fxy_log = log10(fxy);


%Multiply by nchoosek(m+2,2)
lambda_log = sum(fxy_log)

lambda = 10^(lambda_log) ^(nchoosek(m+2,2))

end
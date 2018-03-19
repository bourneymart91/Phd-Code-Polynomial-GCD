function [lambda] = GeometricMean_MatlabMethod(fxy,m,n_k)
% Compute the geometric mean of the entries of the Cauchy like matrix
% C_{n-k}(f(x,y))
%
% % Inputs.
%
% fxy : Coefficients of polynomial f(x,y)
%
% m : Total degree of f(x,y)
%
% n_k : Total degree of v_{k}(x,y)


Cf = BuildSylvesterMatrixPartition_2Polys(fxy, m, n_k);


% Compute geometric mean
lambda = geomean(abs(Cf(Cf~=0)));


end
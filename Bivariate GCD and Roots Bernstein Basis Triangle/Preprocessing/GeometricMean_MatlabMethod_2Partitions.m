function [lambda] = GeometricMean_MatlabMethod_2Partitions(fxy, m, n, o, k)
% Compute the geometric mean of the entries of the Cauchy like matrix
% C_{n-k}(f(x,y)) and C_{o-k}(f(x,y))
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% m : Total degree of f(x,y)
%
% n : (Int) Degree of g(x,y)
%
% o : (Int) Degree of h(x,y)
%
% k : (Int) Index of subresultant matrix
%
% % Outputs
%
% lambda : (Float) Geometric mean of non-zero entries in C_{n-k}(f) and
% C_{o-k}(f)


Cf1 = BuildSylvesterMatrixPartition_2Polys(fxy, m, n - k);
Cf2 = BuildSylvesterMatrixPartition_2Polys(fxy, m, o - k);

% Compute geometric mean
lambda = geomean(...
    [
    abs(Cf1(Cf1~=0)); ...
    abs(Cf2(Cf2~=0)) ...
    ]...
    );


end
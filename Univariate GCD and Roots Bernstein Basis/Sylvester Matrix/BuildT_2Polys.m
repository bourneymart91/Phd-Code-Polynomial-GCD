function T = BuildT_2Polys(fx, gx, k)
%  Build the Toeplitz matrix T = [T1 T2], consisting of coefficients of 
% f(x) and g(x).
%
%
% % Input
%
% fx: (Vector) Coefficients of f(x) in the standard bernstein basis. a_{i}
%
% gx: (Vector) Coefficients of g(x) in the standard Bernstein basis. b_{i}
%
% k : (Int) index of subresultant S_{k} to be formed. (Also degree of GCD)
%
% % Output
%
% T : (Matrix) The partitioned matrix T = [T(f) T(g)].
%

% Get degree of polynomail f
m = GetDegree(fx);
n = GetDegree(gx);

% Build Toeplitz matrix of f, the first partiton.
T1 = BuildT1(fx, n - k);

% Build Toeplitz matrix of g, the second partition.
T2 = BuildT1(gx, m - k);

% Concatenate the partitions.
T = [T1 T2];

end

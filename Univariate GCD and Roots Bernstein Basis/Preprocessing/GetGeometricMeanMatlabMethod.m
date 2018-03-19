function lambda = GetGeometricMeanMatlabMethod(fx,n_k)
% Get geometric mean of the non-zero entries in the matrix C_{n-k}(f(x)),
% which is the (n-k)th order convolution matrix of the polynomial f(x).
% This matrix is also the first partition of S_{k}(f(x),g(x)).
% C_{n-k}(f(x)) may have many forms, "T_{n-k}(f(x))",
% "D^{-1}T_{n-k}(f(x))", "T_{n-k}(f(x))Q_{n-k}" or
% "D^{-1}T_{n-k}(f(x))\hat{Q}" and this is determined by the global
% variable SETTINGS.SYLVESTER_MATRIX_VARIANT
% 
% 
% % Inputs
%
% fx : (Vector) Vector of the coefficients of the polynomial f(x)
%
% n_k : (Integer) Index of the convolution matrix C_{n-k}(f(x))


% Build the partition of the Sylvester matrix
C_fx = BuildSubresultant_Partition_2Polys(fx, n_k);

% Get geometric mean of non-zero entries
lambda = geomean(abs(C_fx(C_fx ~= 0)));

end
function lambda = GetGeometricMeanMatlabMethod_3Polys_3Eqns(fx, n_k, o_k)
% Get the geometric mean of the non-zero entries of the two partitions of
% the kth three-polynomial subresultant matrix, S_{k}(f,g,h), which contain
% the coefficients of the polynomial f(x).
%
% % Inputs
%
% fx : (Vector) The vector of coefficients of the polynomial f(x)
%
% n_k : (Int) n - k : Order of the convolution matrix C_{n - k}(f(x))
%
% o_k : (Int) o - k : Order of the convolution matrix C_{o - k}(f(x))
%
%
% % Outputs
% 
% lambda : (Float) Geometric mean of the non-zero entries of the two
% partitions of the k-th subresultant matrix which contain coefficients of
% the polynomial f(x).


% Build the two partitions of the Sylvester matrix S_{k}(f,g,h) which
% contain the entries of the polynomial f(x).

% Build the first partition C_{n - k}(f(x)) where n is the degree of g(x)
C_f1 = BuildSubresultant_Partition_2Polys(fx, n_k);

% Build the second partition C_{o - k}(f(x)) where o is the degree of h(x)
C_f2 = BuildSubresultant_Partition_2Polys(fx, o_k);

% Get geometric mean of non-zero entries

SetOfEntries = ...
    [...
        abs(C_f1(C_f1~=0))
        abs(C_f2(C_f2~=0))
    ];

lambda = geomean(SetOfEntries);

end
function C1a = BuildDT1Q1_fromPrev(m, n, k, C_0a)
% Build Toeplitz Matrix (Partition of Sylvester Matrix) in the 'Build Up' method.
% m : degree of polynomial f n : degree of polynomial g k : index of
% subresultant S_{k} to be built. C_0a : the preceeding Partition, from
% which C1a will be built.


% Where C_0 is the previous Cauchy Matrix,
k_new = k+1;

% Build Matrix A
A = [zeros((m+n-k),1) diag(1./(1:1:(m+n-k))) ];

% Build Matrix B
Ba = [...
    zeros(1,n-k);...
    diag(1:1:(n-k))
    ];

% Build C1
C1a = A * C_0a * Ba;

%

% if Denominator is included in building toeplitz, then update
% denominator for next S_k
ratio = nchoosek(m+n-(k),n-(k)) ./ nchoosek(m+n-(k_new),n-(k_new));
C1a = C1a .* ratio;


end

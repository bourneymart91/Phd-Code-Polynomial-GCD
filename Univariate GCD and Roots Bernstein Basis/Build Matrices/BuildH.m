function H = BuildH(m,n)
% Build the diagonal matrix H, consisting of two diagonal matrices in block
% diagonal form. H1 and H2. Where H1 consists of binomial coefficients
% corresponding to polynomial f, and H2 consists of binomials corresponding
% to polynomial g. such that H^-1 C(u);C(v) G d = [f;g]
%
%                          Inputs
%
%
% m :   Degree of input polynomial f
%
% n :   Degree of input polynomial g
%
%
% Outputs.
%
% H :   H^{-1} = [H1, 0 ; 0, H2]
% where:
% H1 = diag[1/nchoosek(m,1) ... 1/nchoosek(m,m)]
% H2 = diag[1/nchoosek(n,1) ... 1/nchoosek(n,n)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build H1.
H1 = BuildH1(m);

% Build H2.
H2 = BuildH1(n);

% Get H^{-1} - the block diagonal matrix.
H = blkdiag(H1,H2);

end


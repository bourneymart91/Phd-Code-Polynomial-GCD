function Sk = BuildDTQ_from_prev(S_prev, m, n, k)
% BuildDTQ_from_prev(S_prev, m, n, k)
%
% Construct the kth Sylvester subresultant matrix 
% S_{k}(f,g) = D^{-1}T_{k}(f,g)Q_{k} from the (k-1)th Sylvester 
% subresultant matrix 
%
%
% Inputs.
%
% S_prev : Previous Sylvester matrix S_{k-1}
%
% m : Degree of polynomial f(x,y)
%
% n : Degree of polynomial g(x,y)
% 
% k : Index of Sylvester Matrix To be built.
%
% Outputs.
%
% Sk : The kth Sylvester Subresultant matrix

% Build the matrix VD which premultiplies S_{k-1}(f,g)
VD = BuildVD(m, n, k-1);

% Build the matrix WQ which postmultiplies S_{k-1})(f,g)
WQ = BuildWQ(m, n, k-1);

% Get the Sylvester subresultant matrix S_{k}
Sk = VD * S_prev * WQ;

end


function VD = BuildVD(m, n, k)
% Polynomial f is multiplied by a polynomial g of degree n-k

mat = zeros(m+n-k,m+n-k);

for i1 = 0:1:(m+n-k-1)
    for i2 = 0:1:(m+n-k-1)
        if i1 + i2 <= (m+n-k-1)
            mat(i1+1,i2+1) = 1./((m+n-k)-i1-i2);
        end
    end
end


V1 = (m+n-k) .* mat;

vec_V1 = GetAsVector(V1);

% Get number of zeros
nNonZeros = nchoosek((m+n-k-1)+2,2);

vec_V1 = vec_V1(1:nNonZeros);

V1 = diag(vec_V1);

nRows_ZeroMatrix = nchoosek(m+n-k+1,2);
nCols_ZeroMatrix = nchoosek(m+n-k+1,2)*2/(m+n-k);
zero_mat = zeros( nRows_ZeroMatrix , nCols_ZeroMatrix );

VD = [V1 zero_mat];

end

function WQ = BuildWQ(m, n, k)

part1 = BuildWQ_Partition(n-k);
part2 = BuildWQ_Partition(m-k);

WQ = blkdiag(part1,part2);


end

function WQ = BuildWQ_Partition(n_k)

% Build the WQ matrix for first partition
mat_W1 = zeros(n_k, n_k);

for i1 = 0:1:n_k-1
    for i2 = 0:1:n_k-1
        if i1+i2 <= n_k-1
            mat_W1(i1+1,i2+1) = (n_k-i1-i2) ./ (n_k);
        end
    end
end

vec_W1 = GetAsVector(mat_W1);
nNonZeros_W1 = nchoosek((n_k-1)+2,2);
vec_W1 = vec_W1(1:nNonZeros_W1);
W1 = diag(vec_W1);

nRows_zeroMat = (2/(n_k)) * nchoosek(n_k+1,2);
nCols_zeroMat = nchoosek(n_k+1,2);
WQ = [W1; zeros(nRows_zeroMat, nCols_zeroMat)];

end
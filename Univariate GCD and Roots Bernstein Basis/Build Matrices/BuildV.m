function V = BuildV(m,n,k)
% Build the matrix V
% Used in constructing a Sylvester subresultant matrix S_{k} from S_{k-1}

V = BuildPartitionV(m,n-k);



end

function V1 = BuildPartitionV(m,n_k)
% BuildPartitionV(m,n_k)
%
% Build the matrix V1
% 
% % Inputs
% 
% m : Degree of polynomial f(x)
%
% n_k : Degree of polynomial v(x)


vec = m+n_k+1:-1:1;
vec = (m+n_k+1) ./ vec; 

mat = diag(vec);
mat = [mat zeros(m+n_k+1,1)];

V1 = mat;
    

end


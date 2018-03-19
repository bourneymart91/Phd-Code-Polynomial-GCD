function T1 = BuildT1_Symbolic(m,n_k)
% BuildT1(fxy,m,n_k)
%
% Build the matrix T_{n-k} where T_{n-k}(f(x,y))*v = h. 
%
% The matrix is used in computing the product of f(x,y) and v(x,y).
%
% The matrix T_{n-k}(f) forms part of the convolution matrix 
% C_{n-k}(f) = D^{-1}T_{n-k}(f)Q_{n-k}.
%
% Also used in the Sylvester subresultant matrix S_{k}(f,g) =
% D^{-1}[T_{n-k}(f) T_{m-k}(g)] blkdiag(Q_{n-k} Q_{m-k})
%
% Inputs.
%
% fxy : Coefficients of bivariate polynomial f(x,y) in Bernstein form.
%
% m : Total degree of f(x,y)
%
% n_k : Total degree of v(x,y)
%
% Outputs.
%
% T1 : The partition T_{n-k}(f) of the Sylvester matrix S_{k}(f,g)


fxy = sym('A', [m+1 m+1]);

for i = 0:1:m
    for j = 0:1:m

        if i+j > m
            fxy(i+1,j+1) = 0;
        end
    end
end

disp(fxy)

% Get number of coefficients of v(x,y)
nCoeffs_vxy = nchoosek(n_k+2,2);

% Get number of coefficients in the product fg = h(x,y)
nCoeffs_hxy = nchoosek(m+n_k+2,2);

% Initialise a zero matrix
zero_mat = sym(zeros(m+n_k+1,m+n_k+1));

T1 = sym(zeros(nCoeffs_hxy,nCoeffs_vxy));

% Get fxy with trinomial coefficients
%fxy_tri = GetWithTrinomials(fxy,m);
fxy_tri = fxy;


% Initialise a count
count = 1;

for diag_index = 0:1:n_k
    for i = diag_index:-1:0
        j = diag_index - i;
        
        % Get matrix of coefficients of f(x,y) shifted down by i rows and
        % across by j columns
        temp_mat = zero_mat;
        temp_mat(i+1:i+m+1,j+1:j+m+1) = fxy_tri;
        
        temp_vec = GetAsVector_Symbolic(temp_mat);
        
        % Remove all but the first nchoosek(m+n+2,2) coefficients
        temp_vec = temp_vec(1:nCoeffs_hxy);
        
        % Insert coefficients into the i,j th column of C(f(x,y)).
        
        T1(:,count) = temp_vec;
        
        % Increment counter
        count = count + 1;
        
    end
end



end
function T1 = BuildT1_Total_Bivariate(fxy_matrix, m, n_k)
% BuildT1_Total_bivar(fxy_matrix, m, n-k)
%
% Construct the matrix T(f), which is a partition of the Sylvester matrix
% [T(f) T(g)]. 
% Constructed in terms of total degree. 
% T(f)*v = T(g)*u.
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% m : (Int) Total degree of polynomial f(x,y).
%
% n_k : (Int) Total degree of polynomial v(x,y).
%
% % Outputs
%
% T1 : (Matrix)

% Get size of f(x,y)
[nRows_f, nColumns_f] = size(fxy_matrix);

% Initalise a zero matrix
zero_matrix = zeros(m+n_k+1, m+n_k+1);

% Get number of elements in v(x,y)
nCoefficients_vxy = nchoosek( n_k+1 + 1 , 2);
nColumnsT = nCoefficients_vxy;

% Get the number of elements in f*v
nRowsT = nchoosek((m+n_k+1) + 1,2);

% Initialise a zero matrix.
T1 = zeros(nRowsT, nColumnsT);

% Get number of diagonals in the matrix v(x,y)
nDiagonals_vxy = (n_k+1);

% Initialise a counter
count = 1;

% for each diagonal in M(v)
for tot = 0 : 1 : nDiagonals_vxy - 1
    for i = tot:-1:0
        
        % Initialise a temporary zero matrix
        temp_mat = zero_matrix;
        
        j = tot - i;
        
        if i <= (n_k) && j <= (n_k)
            
            % Wrap the entries of fxy_matrix_padded 1 downward
            % wrap the entires of fxy_matrix 1 place to the right           
            temp_mat((i+1):(nRows_f+i),(j+1):(nColumns_f+j)) = fxy_matrix;
            
            % Produce temporary vector from the coefficients       
            temp_vec = GetAsVector_Version1(temp_mat);
            
            % Remove the last nchoosek(n-k+1,2) entries from temp vec
            % only keep the nchoosek(n-k+2,2) non-zero values which are the
            % coefficients of the polynomial in terms of total degree.
            nCoefficients = nchoosek(m+n_k+2,2);
            temp_vec2 = temp_vec(1:nCoefficients);
            
            % Insert the vector into the matrix T1
            T1(:,count) = temp_vec2;
            
            % Increment counter
            count = count + 1;
        end
    end
    
end
end
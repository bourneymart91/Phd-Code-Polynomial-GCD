function [maximum,minimum] = GetMaxMin2(a_i1i2, i1, i2, m, n_k)
% Note this function assumes the inclusion of Q in the coefficient matrix.
%
% a_{i1i2} : (Float) Coefficient a_{i_{1}, i_{2}} of f(x,y)
%
% i1 : (Int) Index i_{1}
%
% i2 : (Int) Index i_{2}
%
% m : (Int) Degree of f(x,y)
%
% n_k : (Int) Degree of v(x,y)
%
% % Outputs
%
% maximum : (Float) Maximum occurence of the coefficient a_{i_{1},i_{2}} in
% the Sylvester subresultant matrix.
%
% minimum : (Float) Minimum occurence of the coefficient a_{i_{1}, i_{2}}
% in the Sylvester subresultant matrix.


% Build a 2 dimensional vector to store all occurences of the coefficient
% a_{i_{1},i_{2}} from each of the nchoosek(n-k+2,2) columns in S_{k}(f,g).
mat_ai1i2 = zeros(n_k + 1, n_k + 1);


% For each diagonal of v(x,y)
for q = 0 : 1 : n_k
    
    % for each coefficient v_{j1,j2} of v(x,y)
    for j1 = q : -1:0
        
        j2 = q - j1;
        
        entry = GetEntry(a_i1i2, m, n_k, i1, i2, j1, j2);
        
        mat_ai1i2(j1 + 1, j2 + 1) = entry;
        
    end
end

vec = mat_ai1i2(mat_ai1i2~=0);
mat_ai1i2 = vec;


% take absolute values of A
mat_ai1i2 = abs(mat_ai1i2);

try
    [max_row, max_col] = find(mat_ai1i2 == max(mat_ai1i2(:)));
    [min_row, min_col] = find(mat_ai1i2 == min(mat_ai1i2(:)));
    
    % get the maximum and minimum values. Always use (1) since max or min may
    % occur more than once, and we are only interested in one of these values.
    
    maximum = mat_ai1i2(max_row(1), max_col(1));
    minimum = mat_ai1i2(min_row(1), min_col(1));
catch
    maximum = 0;
    minimum = 0;
    
    
end

end



function [entry] = GetEntry(a_i1i2, m, n_k, i1, i2, j1, j2)
%
% % Inputs
%
% a_i1i2 : (Float) Coefficient
%
% m : (Int)
%
% n_k : (Int)
%
% i1 : (Int)
%
% i2 : (Int)
%
% j1 : (Int)
%
% j2 : (Int)

global SETTINGS
switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    
    case 'T'
        entry = a_i1i2 * Trinomial(m, i1, i2);
        
    case 'DT'
        entry = a_i1i2 * Trinomial(m, i1, i2) ./ Trinomial(m + n_k, i1 + j1, i2 + j2);
        
    case 'TQ'
        entry = a_i1i2 * Trinomial(m, i1, i2) * Trinomial(n_k, j1, j2);
        
    case 'DTQ'
        entry = a_i1i2 * Trinomial(m, i1, i2) * Trinomial(n_k, j1, j2) ./ Trinomial(m + n_k, i1 + j1, i2 + j2);
    otherwise
        error('err')
end
end
function [k1k2_pair_mat] = GetPairs_Refined(m,m1,m2,n,n1,n2,k)
% Get possible k1,k2 pairs used in computing the degree of the GCD of two
% polynomials, when the total degree is already given.
%
% % Inputs.
% 
% m : Total degree of polynomial f(x,y)
%
% m1 : Degree of polynomial f(x,y) with respect to x
%
% m2 : Degree of polynomial f(x,y) with respect to y
%
% n : Total degree of polynomial g(x,y) 
%
% n1 : Degree of polynomial g(x,y) with respect to x
%
% n2 : Degree of polynomial g(x,y) with respect to y
%
% k : Degree of d(x,y)

% Initialise the t1,t2 pair matrix
k1k2_pair_mat = [];

% for all possible values of t1
for k1 = k:-1:0
    
    % for all possible values of t2
    for k2 = k:-1:0
        
        % Perform a series of tests
        condA = n1-k1 + n2 -k2 >= n-k;
        condB = m1-k1 + m2 -k2 >= m-k;
        condC = k1 <= n1 && k1 <= m1;
        condD = k2 <= n2 && k2 <= m2;
        condE = n1 - k1 + n2 - k2 <= 2*(n-k);
        condF = m1 - k1 + m2 - k2 <= 2*(m-k);
        condG = n1 - k1 <= n-k;
        condH = n2 - k2 <= n-k;
        condI = m1 - k1 <= m-k;
        condJ = m2 - k2 <= m-k;
        condK = k1+k2 >= k;
        
        
        if condA && condB && condC && condD && condE && condF...
                && condG && condH && condI && condJ && condK
            
            % Get the pair t1 t2 as a new row
            new_row = [k1 k2];
            % Add the new row to the pairs matrix
            k1k2_pair_mat = [k1k2_pair_mat ; new_row];
        end
        
    end
    
    
end


% Remove duplicate rows of the matrix of (t1,t2) pairs
k1k2_pair_mat = unique(k1k2_pair_mat,'rows');



end
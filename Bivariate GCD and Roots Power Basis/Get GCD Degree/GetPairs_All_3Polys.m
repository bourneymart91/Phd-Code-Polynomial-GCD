function k1k2_pair_mat = GetPairs_All_3Polys(m1, m2, n1, n2, o1, o2)
%
% Used in GetGCDDegree
%
% Given that the total degree has been computed, get the set of pairs of
% (k_{1},k_{2}) from which the relative degree is computed.
% Get the set of (k1,k2) pairs for computing S_{k_{1},k_{2}}
%
% % Inputs.
%
% [m1, m2] : Degree of f(x,y) with respect to x and y
%
% [n1, n2] : Degree of g(x,y) with respect to x and y
%
% [o1, o2] : Degree of h(x,y) with respect to x and y
%
% % Outputs.
%
% k1k2_pair_mat : Matrix whose rows consist of pairs k1 and k2.



% The total of t1+t2 must be between t and 2t
k1k2_pair_mat = [];

for k1 = min([m1, n1, o1]):-1:0
    for k2 = min([m2, n2, o2]):-1:0
        
        % Get new pair
        new_row = [k1 k2];
        
        % Add new pair to matrix
        k1k2_pair_mat = [k1k2_pair_mat ; new_row];
        
    end
end



%
% Remove duplicate rows of the matrix of (t1,t2) pairs
k1k2_pair_mat = unique(k1k2_pair_mat,'rows');

%
% if isempty(t1t2_pair_mat)
%
%         condA = n1-k1 + n2 -k2 >= n-t
%         condB = m1-k1 + m2 -k2 >= m-t
%         condC = k1 <= n1 && k1 <= m1
%         condD = k2 <= n2 && k2 <= m2
%         condE = n1 - k1 + n2 - k2 <= 2*(n-t)
%         condF = m1 - k1 + m2 - k2 <= 2*(m-t)
%         condG = n1 - k1 <= n-t
%         condH = n2 - k2 <= n-t
%         condI = m1 - k1 <= m-t
%         condJ = m2 - k2 <= m-t
%         condK = k1+k2 >= t
%         display('')
% end
%


end
function [a,t_degree]=Examples_Univariate_Roots(n)

% This function is a database of pairs of polynomials. The pairs of
% polynomials are indexed by n.

% The matrices a and b define polynomials, where the first column of a
% and b defines the root, and the second column defines the multiplicity
% of the root.

% The matrix c stores the GCD of a and b, in the same format as a and b.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch n
    
    case '1'
        a = [
            0.1     2;
            0.5     1;
            ];
    case '2'
        a = ...
            [
            0.3     1
            ]
        
end


t_degree = sum(a(:,2));

end
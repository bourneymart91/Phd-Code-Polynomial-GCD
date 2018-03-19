%  Sort the terms of a multivariate polynomial f according to 
%  the lexicographical order, and combine the like terms
%
%    Syntax:  >> G = mvPolynSort(F)
%             >> [G,idg] = mvPolynSort(F)
%
%    Input:   F -- (matrix) multivariate polynomial in coeff. matrix
%
%   Output:   G --- (matrix) sorted coeff. matrix of the input polynomial 
%           idg --- (vector) the term indices terms of g
%
%   Example:
%
%      >> F
%      
%      F =
%      
%           0     0     1     3     3     0
%           0     0     1     0     0     1
%           0     0     1     0     0     2
%           8     7    -7     2     4    -5
%      
%
%      >> [g,idg] = mvPolynSort(f)
%
%      G =
%
%           0     3     1     0
%           0     0     1     1
%           0     0     1     2
%          15     6    -7    -5
%
%      idg =
%
%           1     4    14    21

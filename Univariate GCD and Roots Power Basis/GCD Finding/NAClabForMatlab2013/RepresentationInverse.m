%
%  Convert a column vector to the vector in the domain D
%   (the inverse of the function RepresentationVector)
%  
%   Syntax:  v = RepresentationInverse(z, Domain)
%
%    Input:      z -- the column vector representing a vector in the domain
%           Domain -- the vector space domain information represented by
%                      (i) an mxn matrix of 0' and 1's representing
%                          0 entries and variable entries respectively, or
%                     (ii) a polynomial with variable coefficients being 
%                            nonzero, or
%                    (iii) a cell array of matrices in (i) and/or (ii)
%                    assuming the domain is a product of matrix spaces and
%                    polynomial spaces
%

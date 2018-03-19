%
% Find the column vector representation of a vector a in a vector space
%  
%  Syntax:    v = VectorRepresent(a,B)
%
%   Input:    a -- (matrix/string/cell) a (general) vector
%             B -- (cell) basis information, following the output format
%                     Range in HomomorphismMatrix
%             B -- (cell) Basis info of the vector space that is assumed to be
%                     the product of matrix spaces and polynomial spaces, 
%                     each cell entry contains 
%                      (i) an mxn matrix of 0' and 1's representing
%                          zero entries and variable entries respectively,
%                     (ii) or, a 3-dimensional cell with variable names, tuple
%                     degree and term indices for a polynomial component
%  Output:    v -- (vector) column vector representing the vector 'a'
%

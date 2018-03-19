%
%  Solve for the general least squares solution z of L(z) = b 
%     for a linear transformation L
%
%    Syntax:  >> [z, K, lcond, res] = LinearSolve({@L, Domain, para}, b, tol)
%
%    Input:  L  -- Matlab function name for calculating the linear transformation
%                    with syntax  L(X1,...,Xm,Q1,...,Qk) where X1,...,Xm
%                    are variables of L. Each one of X1,...,Xm is 
%                       (i) a matrix (including vector), or
%                      (ii) a polynomial
%                           to be transformed by L.
%                    The remaining parameters Q1,...,Qk must be 
%                        provided as cell entries for para.
%                    If the codomain is a product space, the output of the
%                        function L must be a cell array
%        Domain  -- domain of the linear transformation L represented by
%                      (i) an mxn matrix of 0' and 1's representing
%                          0 entries and variable entries respectively, or
%                     (ii) a polynomial with variable coefficients being 
%                            nonzero, or
%                    (iii) a cell array of matrices in (i) and/or (ii)
%                    assuming the domain is a product of matrix spaces and
%                    polynomial spaces
%          para  -- (cell) parameters needed for running L
%                    **** L must run with >> L(Domain{:},para{:}) ***
%             b  -- (matrix, polynomial or cell) The right-hand side of
%                        the equation L(z) = b, must be in the codomain 
%                        of L
%
%   Output:   z -- the minimum norm solution so that ||L(z)-b|| is minimized
%             K -- an orthonormal basis for the kernel {y | L(y) = 0}
%         lcond -- condition number of the representation matrix.
%           res -- (vector) the residuals of ||L(Z)-b||,||L(K{1}||,...||L(K{n}||
%

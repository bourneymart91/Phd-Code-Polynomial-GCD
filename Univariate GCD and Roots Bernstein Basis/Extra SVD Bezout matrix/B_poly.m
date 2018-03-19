
function [p]=B_poly(A)

% This function implements the convolution operation on the polynomial
% defined by the matrix A, which has two columns, that is, this function
% computes the coefficients of the polynomial defined by A. These are
% the coefficients in the scaled Bernstein basis form of the polynomal.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the number of distinct roots of the polynomial.
r=size(A,1); 
p=1;
 
% Convolve each factor, which is defined by a row of A, separately.
% A(k,1) stores the value of the root, and A(k,2) stores its
% multiplicity.
for k=1:1:r
    w=B_conv(A(k,1),A(k,2));    
    p=conv(p,w);    
end 
    

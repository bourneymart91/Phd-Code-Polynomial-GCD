
function [f_bi]=BuildPoly_Brn(A)
% Obtain polynomial coefficients in the Scaled Bernstein Basis, given a set of roots and multiplicities.
% This function implements the convolution operation on the polynomial
% defined by the matrix A, which has two columns, that is, this function
% computes the coefficients of the polynomial defined by A. These are
% the coefficients in the scaled Bernstein basis form of the polynomial.
% Coefficients of the form ai(m choose i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the number of distinct roots of the polynomial.
r = size(A,1);

% Convolve each factor, which is defined by a row of A, separately.
% A(k,1) stores the value of the root, and A(k,2) stores its multiplicity.

f_bi = 1;
for k = 1:1:r
    w = B_conv(A(k,1),A(k,2));
    f_bi = conv(f_bi,w) ;
end
f_bi = f_bi';
end

function [t]=B_conv(root,mult)
% This function convolves the vector [-r 1-r] with itself m times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs:


% r :   Root

% m :   Multiplicity of root

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Outputs:


% t :   Vector which stores the result from this convolution.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Note that (y-r) = -r(1-y) + (1-r)y and thus the polynomial y-r in the
% power basis is represented as the polynomial -r(1-y) + (1-r)y in the
% scaled Bernstein basis.


if mult==1
    t=[-root,1-root];
else
    
    q=[-root,1-root];
    t=[-root,1-root];
    for k=2:1:mult
        t=conv(t,q);
    end
end
end
function [f]=BuildPolyPwr(A)
% Obtain polynomial coefficients in the power form, given a set of roots 
% and multiplicities.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the number of distinct roots of the polynomial.
r = size(A,1);

% Convolve each factor, which is defined by a row of A, separately.
% A(k,1) stores the value of the root, and A(k,2) stores its multiplicity.

f = 1;
% for each unique root 1,...,r
for k = 1:1:r
    
    w = pwr_conv(A(k,1),A(k,2));
    f = conv(f,w) ;
end
f = f';
end
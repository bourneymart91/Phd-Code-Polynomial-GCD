function [f]=BuildPoly_Pwr(A)
% Obtain polynomial coefficients in the power form, given a set of roots 
% and multiplicities. Coefficients are given in order of descending power.

% The leading coefficient has the highest power. a_{m}x^{m}
% the trailing coefficient is a unit. a_{0}
% Including zero terms
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

% Obtain coefficients so that the leading coefficient has the highest
% power.
f = fliplr(f)

f = f';
end

function [t]=pwr_conv(root,mult)
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
    t=[-root,1];
else
    
    q=[-root,1];
    t=[-root,1];
    for k=2:1:mult
        t=conv(t,q);
    end
end
end

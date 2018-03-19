
function [t]=B_conv(r,m)

% This function convolves the factor vector [-r 1-r] with itself m times.
% The vector t on output stores the result from this convolution.

% Note that (y-r) = -r(1-y) + (1-r)y and thus the polynomial y-r in the
% power basis is represented as the polynomial -r(1-y) + (1-r)y in the 
% scaled Bernstein basis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if m==1
    t=[-r,1-r];
else
    
    q=[-r,1-r];
    t=[-r,1-r];
    for k=2:1:m
        t=conv(t,q);
    end
    
end

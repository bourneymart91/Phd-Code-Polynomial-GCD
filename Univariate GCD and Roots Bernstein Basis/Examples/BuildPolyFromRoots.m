function [f_bi] = BuildPolyFromRoots(root_mult_array)
% B_poly(root_mult_array)
%
% Obtain the coefficients of the polynomial f(x) in scaled Bernstein
% form, from the matrix which contains the roots r_{i} and multiplicities
% m_{i} of each r_{i}.
% 
%
% % Inputs
%
% root_mult_array : Matrix of roots of f(x) and corresponding
%                   multiplicities. where each row consists of a root r_{i}
%                   and its multiplicity m_{i}
%
% % Outputs
%
% f_bi : Coefficients of polynomial f(x) in scaled Bernstein basis.

% Calculate the number of distinct roots of the polynomial.
nRoots = size(root_mult_array,1);

% Convolve each factor, which is defined by a row of A, separately.
% A(k,1) stores the value of the root, and A(k,2) stores its multiplicity.

% Initialise the polynomial
f_bi = 1;

% For each root
for k = 1:1:nRoots
    
    % Get root
    root = root_mult_array(k,1);
    
    % Get multiplicity of the root
    mult = root_mult_array(k,2);
    
    % Get the polynomial of the root convolved with itself m times.
    w = B_conv(root,mult);
    
    % Convolve polynomial f with the polynomial of the root
    f_bi = conv(f_bi,w) ;
    
end

% Transpose to get as a column vector
f_bi = f_bi';


end


function [t] = B_conv(root,mult)
% This function convolves the vector [-r 1-r] with itself m times.
%
% % Inputs
%
% root : root
%
% mult : multiplicity of root
%
% % Outputs
%
% t : vector which stores the result from this convolution.
%


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
function [fx] = BuildPoly_Pwr(root_mult_mat)
% Obtain polynomial coefficients in the power form, given a set of roots 
% and multiplicities. Coefficients are given in order of descending power.
%
% % Inputs
%
% root_mult_mat : (Matrix) Each row contains a root and its multiplicity
%
% % Outputs
%
% fx : (Vector) : Coefficients of polynomial f(x)

% Get the number of distinct factors of the polynomial.
nFactors = size(root_mult_mat,1);

% Convolve each factor, which is defined by a row of A, separately.
% A(k,1) stores the value of the root, and A(k,2) stores its multiplicity.

% Initialise polynomial
fx = 1;

% for each unique root 1,...,r
for k = 1 : 1 : nFactors
    
    % Get temporary polynomial w(x) (x-r)^m
    wx = pwr_conv(root_mult_mat(k,1),root_mult_mat(k,2));
    
    % Multiply f(x) by polynomial w(x)
    fx = conv(fx,wx) ;
    
end

% Obtain coefficients so that the leading coefficient has the highest
% power.
fx = fliplr(fx);

fx = fx';
end

function [poly] = pwr_conv(root, mult)
% This function convolves the vector [-r 1] with itself m times, and
% returns the corresponding polynomial.
% 
%
% % Inputs:
%
% root : (Float) Root
%
% mult : (Int) Multiplicity of root
%
%
% Outputs:
%
% poly : (Vector)  Vector of coefficients of polynomial which is the result 
% from this convolution.
%



% If the multiplicity of the root is only 1, then output the polynomial 
% [ -r , 1]
if mult == 1
    poly = [-root,1];
else
    % Perform polynomial multiplication m times where m is the multiplicity
    % of the root.
    
    % Build polynomial of the simple root
    q = [-root,1];
    
    % Build the starting polynomial 
    poly = [-root,1];
    
    for k = 2 : 1 : mult
        % Multiply the polynomial by the factor. (x-r)
        poly = conv(poly, q);
    end
end

end

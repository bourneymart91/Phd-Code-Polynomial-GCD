
function [vCoefficients]=calculateCoefficients(M,deg,theta)

% This function calculates the coefficients of the GCD of the 
% polynomials f and g from the last non-zero row of the upper triangular 
% matrix M. The degree of the GCD of f and g is deg, and f and g are
% expressed in the Bernstein basis for theta=1, and the modified
% Bernstein basis for theta~=1.

% coe is the vector of coefficients of the GCD in the Bernstein basis,
% and the coefficients are normalised by their 2-norm.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the number of rows of M.
nRows = size(M,1);

% Store the last non-zero row of M in the vector of r and then 
% calculate the coefficients of the GCD.
% The coefficients of the GCD must be scaled by a combinatorial 
% factor.
r = M((nRows - deg),:);
for k = (nRows - deg) : 1 : nRows
    
    r(k)=(by(nRows - 1,k-1)*r(k))/by(deg,(k-nRows+deg));

end

% Extract the coefficients of the GCD from r.
vCoefficients=r(nRows - deg : nRows);

% Transform the coefficients from the modified Bernstein basis to the
% Bernstein basis coefficients.
for k=1:1:length(vCoefficients)
    pk=(2*(nRows-deg))+k-3;
    vCoefficients(k) = vCoefficients(k) / (theta^pk);
end

% Normalise the coefficients by their 2-norm.
normcoeff=norm(vCoefficients);
vCoefficients=vCoefficients/normcoeff;

end
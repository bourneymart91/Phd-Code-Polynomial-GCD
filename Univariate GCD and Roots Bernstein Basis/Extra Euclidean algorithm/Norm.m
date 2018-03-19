
function [p]=Norm(C)

% This function calculates the L2 norm of the coefficients of the
% Bernstein basis polynomial whose coefficients are stored in the 
% vector C. The L2 norm is stored in the scalar p.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the degree of the Bernstein basis polynomial.
n=length(C)-1;

% Calculate the L2 norm of the coefficients of the polynomial.
p=0;

for i=0:1:n
    for j=0:1:n
        p=p+(by(n,i)*by(n,j)*C(i+1)*C(j+1))/(by(2*n,i+j)); 
    end
end

p=p/(2*n+1);
p=sqrt(p);
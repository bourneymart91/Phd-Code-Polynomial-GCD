
function [P]=deg_elevation(p,low_deg,high_deg)

% This function implements the degree elevation operation for the
% Bernstein basis polynomial whose coefficients are stored in the
% vector p.

% p         :  The vector of coefficients of the polynomial that is  
%              to be degree elevated.

% low_deg   :  The degree of the polynomial p.

% high_deg  :  The degree of p after degree elevation.

% P         :  The vector of coefficients of the degree elevated
%              form of p.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=low_deg;  % the degree of p
r=high_deg-n;  % the number of degree elevations

% Initialise the vector of coefficients of the degree elevated
% polynomial.
P=zeros(1,high_deg+1);

% Start the loop for the coefficients of the degree elevated polynomial.
for k=0:1:high_deg
    
    jlow=max(0,k-r);
    jhigh=min(n,k);
    
    for j=jlow:1:jhigh
        t1=(by(r,(k-j))*by(n,j))/by(n+r,k);
        P(k+1)=P(k+1)+(t1*p(j+1));
    end;
    
end;
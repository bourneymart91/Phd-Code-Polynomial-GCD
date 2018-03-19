
function [p]=sub(P)

% This function divides each coefficient of the scaled Bernstein basis 
% polynomial P by its corresponding binomial coefficient in order to 
% obtain the coefficient of the Bernstein basis form p of the polynomial.

% P and p are vectors of the same length.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the degree of the polynomial.
m=length(P)-1;

p=zeros(1,m+1);

% Divide each entry of P by the combinatorial factor in order to
% calculate eachnetry of p.
for k=0:1:m
    p(k+1)=P(k+1)/by(m,k);
end    


function [gm]=GM(c)

% This function calculates the geometric mean gm of the coefficients 
% of the Bernstein basis polynomial whose coefficients are stored in
% the vector c.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

sc=length(c);
c=abs(c);

r=1;
for k=1:1:sc
    r=r*(c(k)^(1/sc));
end

gm=r;


function [B2] = bezout_theta(B,theta)

% This function returns the Bezout matrix B2 for polynomials expressed 
% in the modified Bernstein basis. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% B       :  The Bezout matrix for polynomials expressed in the Bernstein
%            basis.

% theta   :  The value of theta that defines the transformation from
%            the Bernstein basis to the modified Bernstein basis.

% B2      :  The Bezout matrix for polynomials expressed in the 
%            modified Bernstein basis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=size(B,1);
B2=zeros(n,n);

for i=1:1:n
    for j=i:1:n
        B2(i,j)=B(i,j)*(theta^(i+j-2));
        B2(j,i)=B2(i,j);
    end
end


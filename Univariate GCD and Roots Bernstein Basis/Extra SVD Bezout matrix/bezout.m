
function [B] = bezout(p,q)

% This function returns the Bezout matrix B of the Bernstein basis
% polynomials whose coefficients are stored in the vectors p and q. 
% It is assumed that the polynomials are of the same degree.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the degree of p and q, and initialise the Bezout matrix B.
n=length(p)-1;
B=zeros(n,n);

% Construct the Bezout matrix.
for j=1:1:n
    
    for i=j:1:n
        
        if j==1
            B(i,j)=(n/i)*(p(i+1)*q(1)-p(1)*q(i+1));
        elseif (i==n)&&(j~=1)
            B(i,j)=(n/(n-j+1))*(p(i+1)*q(j)-p(j)*q(i+1));
        else    
            B(i,j)=((n*n)/(i*(n-j+1)))*(p(i+1)*q(j)-p(j)*q(i+1))...
                   +(((j-1)*(n-i))/(i*(n-j+1)))*B(i+1,j-1);
        end 
        
        B(j,i)=B(i,j);  % B is symmetric
        
    end
    
end

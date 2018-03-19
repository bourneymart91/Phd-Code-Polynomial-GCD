
function [x] = bs(u)

% This function solves the equation Ax=b by back substitution. The input
% variable u is an upper triangular form of the matrix [A b]. The output 
% variable x is the solution of Ax=b.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n]=size(u);

x=zeros(1,m);
x(m)=u(m,n)/u(m,n-1);

for i=m-1:-1:1
    temp=u(i,n);
    
    for j=m:-1:i+1
        temp=temp-u(i,j)*x(j);
    end
    
    x(i)=temp/u(i,i);    
end
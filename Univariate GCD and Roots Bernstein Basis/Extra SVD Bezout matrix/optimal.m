
function [theta] = optimal(B)

% This function calculates the optimal value of theta for the Bezout
% matrix when the transformation to the modified Bernstein basis is 
% made. The input variable is B, the Bezout matrix for theta=1. The 
% output variable is theta, which is the optimal value of theta.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B=abs(B);
n=size(B,1);

% Calculate the number of constraints in the LP problem whose
% solution is theta.
w=(n+1)*n/2;

% Define the vector in the objective function.
f=[1,-1,0];

Ge=zeros(w,3);
Le=zeros(w,3);
en=zeros(w,1);

k=0;
for i=1:1:n
    for j=i:1:n
        
        k=k+1;
        Ge(k,:)=[1,0,-(i+j-2)];
        Le(k,:)=[0,-1,(i+j-2)];
        en(k)=B(i,j);
            
    end
end

Ge=Ge(1:k,:);
Le=Le(1:k,:);
en=en(1:k,:);


A=(-1)*[Ge;Le];
b=[-log10(en);log10(en)];
x=linprog(f,A,b);

theta=10^x(3);

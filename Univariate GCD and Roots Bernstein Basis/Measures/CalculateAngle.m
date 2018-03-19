
function [a] = CalculateAngle(line,ds)

% This function calculates the angle a between the vector line and the
% space spanned by the columns of the matrix ds.  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transform line to a unit vector.
n=norm(line);
line=line/n;

% Since the angle a is small, it is necessary to use a method that
% is numerically stable. 
[Q,~]=qr(ds,0);
[p,~,~]=svd(line);

[r,c]=size(p);
U2=zeros(r,c-1);
for k=1:1:c-1
    U2(:,k)=p(:,k+1);   
end

v=transpose(U2)*Q;

s=svd(v);
a=s(end);   

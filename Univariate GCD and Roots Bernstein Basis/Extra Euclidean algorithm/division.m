
function [B,r] = division(f,g)

% This function implements the division of two Bernstein basis 
% polynomials. The input vectors f and g store the coefficients of the
% polynomials. The vector B stores the divisor and the vector r
% stores the remainder.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the degree of each polynomial.
deg_f=length(f)-1;
deg_g=length(g)-1;

% Calculate the polynomial of higher degree and assign this polynomial
% to be the dividend (numerator) A. The polynomial of lower degree is
% the divisor (denominator) B.
if deg_f>=deg_g
    m=deg_f;
    A=f;
    n=deg_g;
    B=g;
else
    m=deg_g;
    A=g;
    n=deg_f;
    B=f;
end

e=zeros(m+1,m+2);

% Construct the matrix of linear equations. The matrix Q stores the 
% coefficients of G(x), and the matrix R stores the coefficients of 
% R(x), where F(x)=G(x)Q(x)+R(x).
for k=0:1:m   
    
    Q=zeros(1,m-n+1);
    for j1=max(0,k-n):1:min(m-n,k)
        Q(j1+1)=(by(m-n,j1)*by(n,k-j1)*B(k-j1+1))/by(m,k); 
    end
    
    R=zeros(1,n);
    for j2=max(0,k-m+n-1):1:min(n-1,k)
        R(j2+1)=(by(m-n+1,k-j2)*by(n-1,j2))/by(m,k);
    end
    
    e(k+1,:)=[Q,R,A(k+1)];
    
end

% Calculate the upper triangular form u of e, and then back substitute
% to calculate the vector x that contains the coefficients of the
% quotient q and remainder r.
[~,u]=qr(e);
x=bs(u);
q=x(1:m-n+1);
r=x((m-n+2):m+1);

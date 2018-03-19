
function [resid_u,resid_v] = ScaleCoprime(u,v,f,g,t)

% This function scales the coprime polynomials u and v in order to
% minimise the residual of the equations (u)(d)=f and (v)(d)=g,
% where d is the 
%
%
% Input:
%
% u :- the quotient polynomial, expressed in the Bernstein basis or
%      the modified Bernstein basis.
%
% v :- the quotient polynomial, expressed in the Bernstein basis or
%      the modified Bernstein basis.
%
% f :- the coefficients of f, expressed in the Bernstein basis or
%      the modified Bernstein basis.
%
% g :- the coefficients of g, expressed in the Bernstein basis or
%      the modified Bernstein basis.
%
% t :- the degree of the AGCD.
%
%
% Output
%
%
% resid_u :- the normalised residual of (lambda * u)(d) = f.
%
% resid_x :- the normalised residual of (mu * v)(d) = g.

% The following expression is required

%  {(m-t) choose (i-j)} x {t choose j} 
%  -----------------------------------
%           {m choose i}

% and it is equal to

%  {(m-i) choose (t-j)} x {i choose j} 
%  -----------------------------------
%           {m choose t}

% An identical expression, but with m replaced by n, is also used.



% Construct the right hand side vector = [fx;gx].
bk = [f;g];

% Build the coefficient vector HCG.
C1 = BuildH1C1G(u,t);
C2 = BuildH1C1G(v,t);
        
  
HCG = [C1;C2];

% Obtain the least squares solution by the QR decomposition.
[~,n2] = size(HCG);
[Q,R] = qr(HCG);
R1 = R(1:n2,:);
cd = Q'*bk;
c = cd(1:n2,:);
x_ls = R1\c;
    
dw = x_ls;  % the computed GCD
   
% Calculate the scale factors lambda and mu.
C1dw=C1*dw;
lambda=(C1dw'*f)/(norm(C1dw)^2);

C2dw=C2*dw;
mu=(C2dw'*g)/(norm(C2dw)^2);

% Calculate the residual of each equation.
% The residual of (lambda * u)(d) = f.
num1=(lambda*C1dw)-f;
resid_u=norm(num1)/norm(f);

% The residual of (mu * v)(d) = g.
num2=(mu*C2dw)-g;
resid_v=norm(num2)/norm(g);

end



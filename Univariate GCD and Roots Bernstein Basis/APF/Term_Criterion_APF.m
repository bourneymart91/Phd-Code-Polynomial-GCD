
function [resid_uw,resid_vw] = Term_Criterion_APF(fw,gw,sw,tw,uw,vw,...
    dw,t,alpha)


% This function calculates the termination criterion in the iterative
% solution of the LSE problem when APF is used.

% Input:

% fw,gw     : The polynomials fx and gx. They are expressed in the 
%             Bernstein basis or the modified Bernstein basis,
%             depending on the value of BOOL_PREPROC.

% sw,tw     : The polynomials added to fw and gw.

% uw,vw     : The quotient polynomials in the Bernstein basis or the
%             modified Bernstein basis, depending on the value of
%             BOOL_PREPROC.

% dw        : The GCD, expressed in the Bernstein basis or the
%             modified Bernstein basis.

% t         : The degree of the AGCD.

% alpha     : The value of alpha.

% BOOL_LOG  :  Boolean
%   1 :- Perform combinatorial computations using logarithms.
%   0 :- Perform combinatorial computations using the standard method.

% Output:


% resid_uw  : The normalised residual of the equation (alpha uw)(dw)=fw.
%             The polynomials are expressed in the w variable.

% resid_vw  : The normalised residual of the equation (alpha vw)(dw)=gw.
%             The polynomials are expressed in the w variable.



% Update the polynomials fw and gw.
f = fw + sw;
g = gw + tw;
alphag = alpha*g;

% Construct the matrices H_{1}C_{1}(u)G and H_{2}C_{2}(v)G .
C1 = BuildH1C1G(uw,t);
C2 = BuildH1C1G(vw,t);
          
% Calculate the scale factors lambda and mu.
C1dw=C1*dw;
lambda=(C1dw'*f)/(norm(C1dw)^2);

C2dw=C2*dw;
mu=(C2dw' * alphag)/(norm(C2dw)^2);

% Calculate the residual of each equation.

% (1) The residual of (lambda * uw)(dw) = fw.
num1=(lambda*C1dw)-f;
resid_uw=norm(num1)/norm(f);

% (2) The residual of (mu * vw)(dw) = gw.
num2=(mu*C2dw)-alphag;
resid_vw=norm(num2)/norm(alphag);

end


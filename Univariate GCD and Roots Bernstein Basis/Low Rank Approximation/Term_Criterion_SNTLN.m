
function [resid_ux,resid_vx,resid_uw,resid_vw] = ...
    Term_Criterion_SNTLN(fx,gx,m,n,t,zk,X,alpha,theta,opt_col)


% This function calculates the termination criterion in the iterative
% solution of the LSE problem when SNTLN is used.

% Input:

% fx,gx     : The given inexact Bernstein polynomials that may or may not
%             be normalised by their geometric means, depending in the 
%             value of BOOL_PREPROC. The parameters alpha and theta are
%             not included in fx and gx.

% m,n       : The degrees of m and n.

% t         : The degree of the AGCD.

% zk        : The vector of structured perturbations added to fx and gx.

% X         : The vector of coprime polynomials.

% alpha     : The value of alpha.

% theta     : The value of theta.

% opt_col   : The optimal column.

% BOOL_LOG  :  Boolean
%   1 :- Perform combinatorial computations using logarithms.
%   0 :- Perform combinatorial computations using the standard method.

% Output:

% resid_ux  : The normalised residual of the equation (alpha ux)(dx)=fx.
%             The polynomials are expressed in the x variable.

% resid_vx  : The normalised residual of the equation (alpha vx)(dx)=gx.
%             The polynomials are expressed in the x variable.

% resid_uw  : The normalised residual of the equation (alpha uw)(dw)=fw.
%             The polynomials are expressed in the w variable.

% resid_vw  : The normalised residual of the equation (alpha vw)(dw)=gw.
%             The polynomials are expressed in the w variable.


% Note: This function is only called if BOOL_Q = 1, that is, Q is
% included in the Sylvester matrix.

% Calculate the updated polynomials fx and gx.
fx_output = fx + zk(1:m+1);
gx_output = gx + zk(m+2:end);

% Define the vector vecx from the vector X. This vector
% stores the coprime polynomials.
vecx = [X(1:opt_col-1); -1 ; X(opt_col:end)];
        
% Transform the polynomials fx_output and gx_output to their forms
% in the modified Bernstein basis.
fw_n = GetWithThetas(fx_output,theta);
gw_n = GetWithThetas(gx_output,theta);

% Calculate the coprime polynomials vx and ux. Remove the theta^i
% to transform them to the x variable.
vw = vecx(1:n-t+1);  % the w variable (modified Bernstein basis)
uw = -vecx(n-t+2:end);  % the w variable (modified Bernstein basis)

ux = GetWithoutThetas(uw,theta);
vx = GetWithoutThetas(vw,theta);

% The polynomials fx_output, gx_output, ux and vx are expressed in the
% Bernstein basis. Scale the polynomials ux and vx to minimise the
% errors in the equations (ux)(dx)=fx_output and (vx)(dx)=gx_output. 
[resid_ux,resid_vx] = ScaleCoprime(ux,vx,fx_output,gx_output,t);

% Now repeat the call to ScaleCoprime, but define all the polynomials
% in the modified Bernstein basis. 
[resid_uw,resid_vw] = ScaleCoprime(uw,vw,fw_n,alpha*gw_n,t);

end


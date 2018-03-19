
function [r] = CalculateResidualSVD(ck,Ak)
%% Calculate Residual by Psuedo Inverse Method
% This function calculates the residual r of an approximate linear
% algebraic equation whose coefficient matrix is Ak and right hand 
% side vector is ck. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain x, the approximate solution vector 
x = pinv(Ak)*ck;

% Calculate residual
r = ck-(Ak*x);

% Normalise the residual
r = norm(r);

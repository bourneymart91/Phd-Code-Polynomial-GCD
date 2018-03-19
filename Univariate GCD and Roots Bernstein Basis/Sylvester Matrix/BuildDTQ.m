function [DTQ] = BuildDTQ(fx, gx, k)
% BuildDTQ(fx,gx,t)
%
% Build the matrix DTQ = D^{-1}T(f,g)*Q.
%
% Inputs.
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% gx : (Vector) Coefficients of polynomial g(x)
%
% k : (Int) index of subresultant
%
% % Outputs
%
% DTQ 

% Get degree of polynomial f(x)
m = GetDegree(fx);

% Get degree of polynomial g(x)
n = GetDegree(gx);

% Build matrix D^{-1}
D = BuildD_2Polys(m, n - k);

% Build matrix T(f,g) = T_{n-k}(f) T_{m-k}*(g)
T = BuildT_2Polys(fx, gx, k);

% Build matrix Q = [Q_{n-k} Q_{m-k}]
Q = BuildQ_2Polys(m, n, k);

% Get D^{-1} * T(f,g) * Q
DTQ = D*T*Q;



end
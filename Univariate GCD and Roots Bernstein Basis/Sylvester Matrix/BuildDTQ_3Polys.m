function [DTQ] = BuildDTQ_3Polys(fx, gx, hx, k)
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
% hx : (Vector) Coefficients of polynomial h(x)
%
% k : (Int) index of subresultant
%
% % Outputs
%
% DTQ 

% Get degree of polynomial f(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

% Build matrix D^{-1}
D = BuildD_3Polys(m, n - k,o - k);

% Build matrix T(f,g) = T_{n-k}(f) T_{m-k}*(g)
T = BuildT_3Polys(fx, gx, hx, k);

% Build matrix Q = [Q_{n-k} Q_{m-k}]
Q = BuildQ_3Polys(m, n, o, k);

% Get D^{-1} * T(f,g) * Q
DTQ = D*T*Q;



end
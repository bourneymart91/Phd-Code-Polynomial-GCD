function T = BuildT_3Polys_3Eqns(fx, gx, hx, k)
%  Build the Toeplitz matrix T = [T1 T2], consisting of coefficients of 
% f(x) and g(x).
%
%
% % Input
%
% fx : (Vector)
%
% gx : (Vector)
%
% hx : (Vector)
%
% k : (Int)

% Get the degree of polynomail f(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

% % Build the block diagonal matrix of T_{n-k}(f) and T_{o-k}(f)
T1 = BuildT1(fx, n - k);
T2 = zeros(m + n - k + 1, o - k + 1);
T3 = BuildT1(gx, m - k);
T4 = zeros(m + o - k + 1, n - k + 1);
T5 = BuildT1(fx, o - k);
T6 = BuildT1(hx, m - k);
T7 = BuildT1(hx, n - k);
T8 = BuildT1(gx, o - k);
T9 = zeros(n + o - k + 1, m - k + 1);





% Concatenate the partitions.
T = [T1 T2 T3 ; T4 T5 T6; T7 -T8 T9];


end

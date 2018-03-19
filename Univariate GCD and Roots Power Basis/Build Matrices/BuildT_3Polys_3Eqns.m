function Sk = BuildT_3Polys_3Eqns(fx, gx, hx, k)
% Build the Sylvester Subresultant matrix S_{k}(f,g)
%
% % Inputs
%
% fx : Coefficients of the polynomial f(x)
%
% gx : Coefficients of the polynomial g(x)
%
% hx : Coefficients of the polynomial h(x)
%
% k : Index of Sylvester Subresultant matrix to be constructed.


% Get degree of polynomial f(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

T1 = BuildT1(fx, n - k);
T2 = zeros( m + n - k + 1, o - k + 1); 
T3 = BuildT1(gx, m - k);

T4 = zeros( m + o - k + 1, n - k + 1);
T5 = BuildT1(fx, o - k);
T6 = BuildT1(hx, m - k);

T7 = BuildT1(hx, n - k);
T8 = BuildT1(gx, o - k);
T9 = zeros( n + o - k + 1, m - k + 1);

Sk = [ T1 T2 T3; T4 T5 T6; T7 -T8 T9];

end
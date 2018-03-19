function dx = EuclideanAlgorithm(fx,gx)
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x)
%
% % Outputs
%
% dx : (Vector) Coefficients of polynomial d(x)
 

% Convert the polynomials to the power basis
fx_pwr = BernsteinToPower(fx);
gx_pwr = BernsteinToPower(gx);

a = fx_pwr;
b = gx_pwr;

q = 0
r = a;
d = GetDegree(b);
c = b(1);

while GetDegree(r) >= d
    s = r(1)



end
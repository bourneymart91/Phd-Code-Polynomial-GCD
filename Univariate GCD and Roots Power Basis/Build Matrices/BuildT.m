function Sk = BuildT(fx, gx, k)
% Build the Sylvester Subresultant matrix S_{k}(f,g)
%
% % Inputs
%
% fx : Coefficients of the polynomial f(x)
%
% gx : Coefficients of the polynomial g(x)
%
% k : Index of Sylvester Subresultant matrix to be constructed.


% Get degree of polynomial f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Build the matrix T_{n-k}(f(x)) and buiild the matrix T_{m-k}(g(x))
C1 = BuildT1(fx, n - k);
C2 = BuildT1(gx, m - k);

% Build the partitioned matrix S_{k}(f(x),g(x))
Sk = [C1 C2];

end
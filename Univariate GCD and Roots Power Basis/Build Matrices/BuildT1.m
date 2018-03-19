function [T_f] = BuildT1(fx,n_k)
% Given the polynomial f(x) build the partition C1 of the Sylvester matrix.
% C(f) * v = C(g)*u.
%
% Inputs
%
% fx : Coefficients of the polynomial f(x)
%
% n_k : Degree of the polynomial v(x)
%

% Get degree of polynomial f(x)
m = GetDegree(fx);

% Initialise an empty matrix T_{n-k}(f(x))
T_f = zeros(m + n_k + 1, n_k + 1);

% Build the matrix T_{n-k}(f(x))
for i = 0:1:n_k
   
    T_f(i + 1 : i + m + 1, i + 1) = fx;
   
end

end
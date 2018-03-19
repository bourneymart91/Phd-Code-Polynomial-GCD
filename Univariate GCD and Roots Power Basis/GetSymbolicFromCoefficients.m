function fx_sym = GetSymbolicFromCoefficients(fx)
% Given the vector of coefficients of a polynomial, return a symbolic
% expression representing the polynomial
% 
% Inputs.
%
% fx : Vector of coefficients of polynomial f(x)
% 
% Outputs.
%
% fx_sym : Symbolic expression representing polynomial
%
% eg. GetSymbolicFromCoefficients([1;2;3;4])
%   --> '1 + 2*x + 3*x^2 + 4*x^3'


% Intialise the symbolic variable x
x = sym('x');

% Get the degree of f(x)
m = GetDegree(fx);

% Get symbolic vector
sym_vector = diag(x.^(0:1:m)) * fx;

% Get the symbolic expression
fx_sym = sum(sym_vector);

end




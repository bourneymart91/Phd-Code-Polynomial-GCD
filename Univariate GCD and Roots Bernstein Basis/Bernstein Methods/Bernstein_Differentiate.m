function gx = Bernstein_Differentiate(fx)
% Differentiate the polynomial f(x) whose coefficients are given in the
% Bernstein basis. Let the derivative of f(x) be known as g(x).
%
% Inputs.
%
% fx : Coefficients of input polynomial f(x)
%
% Outputs.
%
% gx : Deriveative of polynomial f(x).

% Get degree of polynomial f(x)
m = GetDegree(fx);

% Initialise the vector of the derivative of f
gx = zeros(m, 1);


%Loop through all coefficients of g(x)
for i = 0:1:m - 1
    gx(i+1) = m * (fx(i+1+1) - fx(i+1));
end




end
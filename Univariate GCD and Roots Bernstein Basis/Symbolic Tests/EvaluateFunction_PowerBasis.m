function [f] = EvaluateFunction_PowerBasis(a,b,inc,fx)
%% Evaluate the function fx over interval [a,b] where coefficients fx are defined in the power basis.
%
%
%                            Inputs.
%
%
% a - lower limit of interval
%
% b - upper limit of interval
%
% inc - size of steps within the interval [a,b]
%
% fx - Coefficients of Polynomial f in the power basis.
%
% Outputs.
%
% f - vector of evaluated points of f on interval [a:inc:b]





% Get set of evaluation points, evenly distributed over interval [a,b]
x = a:inc:b;

% Get set of y values where y = f(x).
f = polyval(fx,x);


end

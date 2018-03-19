function [ft] = Bernstein_Evaluate(fx,t)
% Given the function f(t), evaluate f(t) a point t.
%
% Inputs
%   
%   f   :   Coefficients of polynomial f(t) in Bernstein form.
%
%   t   :   Point at which to evaluate f(t)


% Get degree of polynomial f(x).
m = size(fx,1) - 1;

% Initialise the evaluation sum.
sum = 0;

% For each coefficient of f(x) evaluate.
for i = 0:1:m
    
    sum = sum + fx(i+1) .* nchoosek(m,i) .* ((1-t).^(m-i)) .* (t^i);
    
end

ft = sum;

end



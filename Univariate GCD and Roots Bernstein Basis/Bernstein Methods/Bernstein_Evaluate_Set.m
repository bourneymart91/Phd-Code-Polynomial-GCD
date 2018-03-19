function [f] = Bernstein_Evaluate_Set(a,b,inc,fx)
% Evaluate the function f(x) over interval [a,b].
%
%
%                           Inputs.
%
% a :   lower limit of interval
%
% b :   upper limit of interval
%
% inc : size of steps within the interval [a,b]
%
% fx :  Coefficients of Polynomial f in the bernstein basis.
%
% f :   Vector of values corresponding to the evaluated function f, at 
% points x = [a:inc:b]
%
%

i = 0;


x = a:inc:b;

f = zeros(1,length(x));

for c = a:inc:b
    i = i+1;
    f(i) = Bernstein_Evaluate(fx,c);
end

end

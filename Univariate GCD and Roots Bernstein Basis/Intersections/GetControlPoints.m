function Pk = GetControlPoints(a,b,f)
% Get set of control points Pk for the polynomial f in Bernstein form over
% the interval [a,b]
%
% Inputs.
%
%
% a :- Interval Lower Limit.
%
% b :- Interval Upper Limit.
%
% f :- Coefficients of Polynomial in Bernstein form
%
% Outputs.
%
%
% Pk :- Set of control points.

%%
% Get degree of polynomial f
m = GetDegree(f);

% Initialise the matrix of control points.
% Number of control points = number of coefficients
% 2 columns, m+1 rows, where m = degree of polynomial f(x)

Pk = [m+1,2];

% for each control point, assign value.
for i = 0:1:m
    Pk(i+1,:) = [a+(i/m).*(b-a)    f(i+1)];
end



end
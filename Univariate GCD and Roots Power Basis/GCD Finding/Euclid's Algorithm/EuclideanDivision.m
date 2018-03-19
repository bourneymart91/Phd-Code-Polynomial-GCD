function [] = EuclideanDivision(a,b)
%
%
% Inputs
%
% a : Polynomial in the variable x
%
% b : Polynomial in the variable x
%
% Outputs
%
% q : the quotient
%
% r : remainder

q = 0;
r = a;
d = deg(b);
c = lc(b);
while deg(r) >= d

    s = lc(r) ./ c 
    
end

end

function m = deg(a)

[r,~] = size(a);
m = r - 1;

end

function coef = lc(a)
% The leading coefficient is the last coefficient in the vector.

coef = a(end);

end
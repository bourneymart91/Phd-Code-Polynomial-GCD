function [dist] = GetDistanceBetweenPolynomials(fxy_exact,fxy_computed,name)
% Given the two polynomials f(x) exact and f(x) computed. Calculate the
% distance between them.
%
% % Inputs
%
% fxy_exact : Matrix of coefficients of a polynomial in its exact form
%
% fxy_computed : Matrix of coefficients of a polynomial in its computed
% form.
%
% name : Name of the polynomial eg. 'd(x,y)'

% Ensure that f(x,y) exact and f(x,y) computed are both in term of
% total degree.
fxy_computed = Normalise(fxy_computed);
fxy_exact = Normalise(fxy_exact);



%fxy_computed - fxy_exact

try
    dist = norm(fxy_exact - fxy_computed) ./ norm(fxy_exact);
catch
    dist = 10000;
end


fprintf([mfilename ' : ' sprintf('Distance of %s exact from %s computed a - b / a : \t %2.4e \n',name,name,dist)]);







end

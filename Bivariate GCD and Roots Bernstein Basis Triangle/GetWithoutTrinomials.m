function fxy = GetWithoutTrinomials(fxy_tri, m)
% Given the coefficients of the polynomial f(x,y) in Bernstein form, get
% the coefficients in the scaled Bernstein form. That is, coefficients of
% f(x,y) with trinomials included.
%
% % Inputs
%
% fxy_tri : Coefficients of polynomial f(x,y)
%
% m : Degree of polynomial f(x,y)
%
% % Outputs
%
% fxy_tri : Coefficients of polynomial f(x,y) in scaled Bernstein basis.


% Initialise a zero matrix to store coefficients of f(x,y)
fxy = zeros(m + 1, m + 1);

for i = 0 : 1 : m
    for j = 0 : 1 : m - i
        
        fxy(i + 1, j + 1) = fxy_tri(i + 1, j + 1) ./ Trinomial(m, i, j);
        
    end
end

end
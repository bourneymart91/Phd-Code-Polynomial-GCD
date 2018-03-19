function [fxy_bb] = PowerToBernstein(fxy, m)
%
% Given the set of coefficients a_{i,j} of polynomial f(x,y), return the 
% coefficients b_{i,j} of the same polynomial in Bernstein form 
%
% Power form.
% a_{i,j} x^{i}y^{j}
%
% Bernstein form.
% b_{i,j} B_{i,j}(x,y)
%
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y). 
%
% m : (Int) Total degree of f(x,y)
%
% % Outputs
%
% fxy_bb : (Matrix) Coefficients of polynomial f(x,y) in the Bernstein
% basis.

% Get the degree of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

temp = zeros(m+1, m+1);
temp(1 : m1+1, 1 : m2+1) = fxy;
fxy = temp;

if m == 0
   fxy_bb = 1;
   return
end

b = zeros(m+1, m+1);

% i and j denote the row and column index for the coefficient matrix of the
% polynomial in bernstein form.
for i = 0:1:m
    for j = 0:1:m-i
        
        b_ij = 0;
        
        for t = 0:1:j
            for s = 0:1:i
                b_ij = b_ij + ...
                    (...
                        nchoosek(i,s) ...
                        * nchoosek(j,t) ...
                        / Trinomial(m,s,t)...
                     )...
                     *fxy(s+1,t+1);
            end
        end
        
        b(i+1,j+1) = b_ij;
    end
end

fxy_bb = b;

end
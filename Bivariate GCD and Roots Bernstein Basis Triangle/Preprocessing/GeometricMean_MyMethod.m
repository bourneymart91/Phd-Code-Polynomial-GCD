function [lambda] = GeometricMean_MyMethod(fxy, m, n_k)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) in scaled Bernstein form
%
% m : (Int) Total degree of f(x,y)
%
% n_k : (Int) Total degree of polynomial v(x,y)
%
% % Outputs
%
% lambda : (Float) Geometric mean of entries of the matrix C_{n-k}(f)


myprod = 1;
count = 0;

fxy = abs(fxy);


for q = 0:1:m
    for i1 = 0:1:q
        for p = 0:1:n_k
            for j1 = 0:1:p
   
                i2 = q - i1;
                j2 = p - j1;
                
                ai1i2 = fxy(i1+1,i2+1);
                
                myprod = myprod .* ((ai1i2 * Trinomial(m,i1,i2) * Trinomial(n_k,j1,j2)) ./ Trinomial(m+n_k,i1+j1,i2+j2));
                count = count + 1;
            end
        end
    end
end

nEntries_fxy = nchoosek(m+2,2);

nNonZeroEntries_Cf = nEntries_fxy * nchoosek(n_k+2,2);

lambda = myprod .^(1./(nNonZeroEntries_Cf));


end
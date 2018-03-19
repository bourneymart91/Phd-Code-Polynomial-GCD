function [dist] = GetMatrixDistance(fxy_exact, gxy_approx)
%
% % Inputs
%
% fxy : (Matrix) Exact Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Approximated Coefficients of polynomial g(x,y)
%
%
% % Outputs
%
% dist : Distance between two polynomials

fxy_exact = fxy_exact./norm(fxy_exact);
gxy_approx = gxy_approx./norm(gxy_approx);


fxy_exact = MakePositive(fxy_exact);
gxy_approx = MakePositive(gxy_approx);

try
    
    dist = norm(fxy_exact - gxy_approx) ./ norm(fxy_exact);
    
catch
    dist = 1000;
end

end


function fxy = MakePositive(fxy)
% Maket the leading coefficient positive


if(sum(sum(fxy)) < 0)
   
    fxy = fxy .* -1;
    
end

end
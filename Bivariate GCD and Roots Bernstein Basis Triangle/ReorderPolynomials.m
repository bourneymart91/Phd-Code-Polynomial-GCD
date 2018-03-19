function [fxy, gxy, hxy, uxy, vxy, wxy, dxy, m, m1, m2, n, n1, n2, o, o1, o2,...
    t, t1, t2] = ReorderPolynomials(fxy, gxy, hxy,...
    uxy, vxy, wxy, dxy, m, m1, m2, n, n1, n2, o, o1, o2, t, t1, t2)
%
% % Inputs
%
%


% Get polynomial with minimum degree
[~, index] = min([m n o]);


switch index
    
    case 2 % Polynomial g has minimum degree
        
        [fxy, gxy] = SwapPolys(fxy, gxy);
        [uxy, vxy] = SwapPolys(uxy, vxy);
        [m,n] = SwapDegree(m,n);
        
        
    case 3 % Polynomial h has minimum degree
        
        [fxy, hxy] = SwapPolys(fxy, hxy);
        [uxy, wxy] = SwapPolys(uxy, wxy);
        [m,o] = SwapDegree(m,o);
        
    otherwise
        % Do nothing
end


end

function [fxy, gxy ] = SwapPolys(fxy, gxy)


fxy_temp = gxy;
gxy = fxy;
fxy = fxy_temp;


end

function [m, n] = SwapDegree(m, n)


m_temp = n;
n = m;
m = m_temp;


end


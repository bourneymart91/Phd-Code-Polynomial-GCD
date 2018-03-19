function Sk = BuildSubresultant_3Polys(fx, gx, hx, k)
% Build the Sylvester Subresultant matrix S_{k}(f,g)
%
% % Inputs
%
% fx : Coefficients of the polynomial f(x)
%
% gx : Coefficients of the polynomial g(x)
%
% hx : Coefficients of the polynomial h(x)
%
% k : Index of Sylvester Subresultant matrix to be constructed.


% Get degree of polynomial f(x), g(x) and h(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

global SETTINGS

switch SETTINGS.N_EQUATIONS_SYLVESTER_MATRIX
    
    case '2'
        Sk = BuildT_3Polys_2Eqns(fx, gx, hx, k);
        
        
        
    case '3'
        Sk = BuildT_3Polys_3Eqns(fx, gx, hx, k);
        
    otherwise
        error('err')
end
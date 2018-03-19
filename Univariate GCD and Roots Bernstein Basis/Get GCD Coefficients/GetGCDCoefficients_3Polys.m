function [dx] = GetGCDCoefficients_3Polys(ux, vx, wx, fx, gx, hx, k)
% Get The Coefficients of the approximate GCD using Quotient Polynomials.
%
% % Inputs
%
% ux : (Vector) The coefficients of the polynomial u(x)
%
% vx : (Vector) The coefficients of the polynomial v(x)
%
% wx : (Vector) The coefficients of the polynomial w(x)
%
% fx : (Vector) The coefficients of the polynomial f(x)
%
% gx : (Vector) The coefficients of the polynomial g(x)
%
% hx : (Vector) The coefficients of the polynomial h(x)
%
% k : (Int) The degree of common divisor.
%
%
% % Outputs
%
%
% dx : (Vector) The coefficients of the approximation of the polynomial d(x)

if (nargin ~= 7)
   error('Not enough input arguments'); 
end

% Global variables
global SETTINGS



switch SETTINGS.GCD_COEFFICIENT_METHOD
    
    case 'ux and vx'
        % Build solution vector bk = [f; g; h]
        bk = [fx ; gx; hx];
        
        % Build the coefficient vector HCG
        HCG = BuildHCG_3Polys(ux, vx, wx, k);
        
        % Get the vector d(w), which is the solution of a problem of the form Ax=b
        dx = SolveAx_b(HCG, bk);
        
        
    case 'ux'
        bk = fx;
        
        H1C1G = BuildH1C1G(ux, k);
        
        dx = SolveAx_b(H1C1G, bk);
        
        
    otherwise
        error('GCD_COEFFICIENT_METHOD is either (ux) or (ux and vx)')
end



end


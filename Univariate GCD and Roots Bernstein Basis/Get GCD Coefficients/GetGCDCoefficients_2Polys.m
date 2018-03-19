function [dx] = GetGCDCoefficients_2Polys(ux, vx, fx, gx, t)
% Get The Coefficients of the approximate GCD using Quotient Polynomials.
%
% % Inputs
%
% ux : (Vector) Coefficients of the polynomial u(x) 
%
% vx : (Vector) Coefficients of the polynomial v(x) 
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x) 
%
% k : (Int) Degree of common divisor of f(x) and g(x)
%
% alpha : (Float) Optimal value of \alpha
%
% theta : (Float) Optimal value of \theta


% Global variables
global SETTINGS

% Coefficients of d(x) are computed from the system : 
%
% [C(u)] * d = f 
%
% or 
%
% [C(u) ; C(v) ] * d = [f ; g]

switch SETTINGS.GCD_COEFFICIENT_METHOD
    
    case 'ux and vx'
        
        % Build solution vector bk = [f;g]
        bk = [fx ; gx];
        
        % Build the coefficient vector HCG
        HCG = BuildHCG_2Polys(ux, vx, t);
        
        % Get the vector d(w), which is the solution of a problem of the form Ax=b
        dx = SolveAx_b(HCG,bk);
        
    case 'ux'
        
        bk = fx;
        
        H1C1G = BuildH1C1G(ux, t);
        
        dx = SolveAx_b(H1C1G, bk);
        
        
    otherwise
        error('GCD_COEFFICIENT_METHOD is either (ux) or (ux and vx)')
end





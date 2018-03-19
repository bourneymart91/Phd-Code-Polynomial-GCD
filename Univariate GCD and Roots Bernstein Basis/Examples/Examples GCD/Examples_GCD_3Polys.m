function [fx, gx, hx, dx, ux, vx, wx] = Examples_GCD_3Polys(ex_num, ex_num_variant)
% Inputs. 
%
% ex_num : (String) Example number
%
% ex_num_variant : (String) 'a', 'b' or 'c' determines the ordering of the
%   polynomials found in the example file.
%
% Outputs.
%
% fx : (Vector) The coefficients of the polynomial f(x)
%
% gx : (Vector) The coefficients of the polynomial g(x)
%
% hx : (Vector) The coefficients of the polynomial h(x)
%
% dx : (Vector) The coefficients of the polynomial d(x), the GCD of f(x) 
%       g(x) and h(x)
%
% ux : (Vector) The coefficients of the polynomial u(x), given by f(x)/d(x)
%
% vx : (Vector) The coefficients of the polynomial u(x), given by g(x)/d(x)
%
% wx : (Vector) The coefficients of the polynomial u(x), given by h(x)/d(x)



% Get the coefficients of the polynomials f(x), g(x), h(x) and the GCD 
% d(x), as well as the coefficients of the cofactor polynomials u(x), v(x) 
% and w(x).
[fx, gx, hx, dx, ux, vx, wx] = ...
    Examples_GCD_FromCoefficients_3Polys(ex_num, ex_num_variant);
        
        


end




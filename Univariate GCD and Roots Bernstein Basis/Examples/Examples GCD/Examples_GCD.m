function [fx, gx, dx, ux, vx] = Examples_GCD(ex_num)
% Inputs. 
%
% ex_num : (String) Example number
%
% Outputs.
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x)
%
% dx : (Vector) Coefficients of the polynomial d(x), the GCD of f(x) and g(x)
%
% ux : (Vector) Coefficients of the polynomial u(x), given by f(x)/d(x)
%
% vx : (Vector) Coefficients of the polynomial v(x), given by g(x)/d(x)


[fx, gx, dx, ux, vx] = Examples_GCD_FromCoefficients(ex_num);
        
    



end




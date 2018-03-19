function dx = GetGCDCoefficients_3Polys(ux, vx, wx, fx, gx, hx, k)
% Get the common divisor of degree k, of the polynomials f(x) and g(x).
%
% % Inputs
%
% ux : Coefficients of the polynomial u(x)
%
% vx : Coefficients of the polynomial v(x)
%
% wx : Coefficients of the polynomial w(x)
%
% fx : Coefficients of the polynomial f(x)
%
% gx : Coefficients of the polynomial g(x)
%
% hx : Coefficients of the polynomial h(x)
%
% k : Degree of the common divisor d(x)
%
% % Outputs
%
% dx : Coefficients of the polynomial d(x) of degree k.


global SETTINGS

switch SETTINGS.GCD_COEFFICIENT_METHOD 
    case 'ux'
        
        dx =  GetGCD_ux(ux,fx,k);
        
    case 'ux and vx'
        
        dx = GetGCD_ux_and_vx(ux, vx, wx, fx, gx, hx, k);
        
end


end



function [dx] = GetGCD_ux(ux,fx,k)
%
% % Inputs
%
% ux : Coefficients of the polynomial u(x)
%
% vx : Coefficients of the polynomial v(x)
%
% k : Degree of the common divisor

% Build the toeplitz matrix C_{k}(u)
C_u = BuildT1(ux,k);

% Build the Right hand side vector [f(\omega) ; alpha.* g(\omega)]
rhs_vec = fx;

% Get the vector d.
d_ls = SolveAx_b(C_u,rhs_vec);

% Calculate d(\omega)
dx = d_ls;

end


function [dx] = GetGCD_ux_and_vx(ux, vx, wx, fx, gx, hx, k)
%
% % Inputs
%
% ux : Coefficients of polynomial u(x)
%
% vx : Coefficients of polynomial v(x)
%
% wx : Coefficients of polynomial w(x)
%
% fx : Coefficients of polynomial f(x)
%
% gx : Coefficients of polynomial g(x)
%
% hx : Coefficients of polynomial h(x)
%
% k : Degree of common divisor d(x)

% Get the GCD d(x) by the APF
C_u = BuildT1(ux,k);
C_v = BuildT1(vx,k);
C_w = BuildT1(wx,k);


C1 = ...
    [
    C_u;
    C_v;
    C_w;
    ];

% Build the Right hand side vector [f(\omega) ; alpha.* g(\omega)]
rhs_vec = ...
    [
    fx;
    gx;
    hx;
    ];

% Get the vector d.
d_ls = SolveAx_b(C1,rhs_vec);

% Calculate d(\omega)
dx = d_ls;


end
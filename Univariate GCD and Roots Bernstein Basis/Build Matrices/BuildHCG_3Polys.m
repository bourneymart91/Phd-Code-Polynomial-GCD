function [HCG, H1C1G, H2C2G,H3C3G] = BuildHCG_3Polys(ux, vx, wx, t)
% Build the matrix HCG, such that H*C(u,v)*G * d = [f;g]
%
% % Inputs
%
% [ux, vx, wx] : Coefficients of the polynomial u(x), v(x) and w(x)
%
% t : Degree of GCD d(x,y)
%
% % Outputs
%
% HCG :
%
% H1C1G :
%
% H2C2G :
%
% H3C3G :



% Build the first part
H1C1G = BuildH1C1G(ux,t);

H2C2G = BuildH1C1G(vx,t);

H3C3G = BuildH1C1G(wx,t);

HCG = ...
    [
        H1C1G;
        H2C2G;
        H3C3G;
    ];

end
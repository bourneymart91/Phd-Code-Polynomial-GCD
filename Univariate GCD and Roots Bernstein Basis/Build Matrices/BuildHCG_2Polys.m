function [HCG, H1C1G, H2C2G] = BuildHCG_2Polys(ux, vx, t)
% Build the matrix HCG, such that H*C(u,v)*G * d = [f;g]
%
% % Inputs
%
% [ux, vx] : Coefficients of the polynomial u(x) and v(x)
%
% t : Degree of GCD d(x,y)
%
%
% % Outputs
%
% HCG :
%
% H1C1G :
% 
% H2C2G :


% Build the first part
H1C1G = BuildH1C1G(ux,t);

H2C2G = BuildH1C1G(vx,t);

HCG = ...
    [
        H1C1G;
        H2C2G
    ];

end
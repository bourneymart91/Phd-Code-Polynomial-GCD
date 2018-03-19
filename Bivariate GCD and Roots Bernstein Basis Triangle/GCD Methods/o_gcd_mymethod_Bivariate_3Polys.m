function [fxy_o, gxy_o, hxy_o, dxy_o, uxy_o, vxy_o, wxy_o, t, rank_range] = ...
    o_gcd_mymethod_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, limits_t, rank_range)
% o_gcd_mymethod(fxy, gxy, m, n, limits_t)

% Given two bivariate polynomials, return the GCD d(x,y) and the coprime
% polynomials u(x,y) and v(x,y) where
%
%
% f/d = u
% g/d = v
%
% fv-gu = 0
%
%
% % Inputs.
%
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of the polynomial h(x,y)
%
% m : (Int) Total degree of the polynomial f(x,y)
%
% n : (Int) Total degree of the polynomial g(x,y)
%
% o : (Int) Total degree of the polynomial h(x,y)
%
%
% % Outputs
%
% uxy : (Matrix) Coefficients of the polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of the polynomial v(x,y)
%
% wxy : (Matrix) Coefficients of the polynomial w(x,y)
%
% dxy : Calculated coefficients of the GCD d(x,y)

% Compute the degree of the GCD
[t, GM_fxy, GM_gxy, GM_hxy, alpha, beta, gamma, th1, th2, rank_range] = ...
    GetGCDDegree_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, limits_t, rank_range);

LineBreakMedium();
fprintf([mfilename sprintf('Geometric Mean f(x,y) 1 : %f \n', GM_fxy)]);
fprintf([mfilename sprintf('Geometric Mean g(x,y)   : %f \n', GM_gxy)]);
fprintf([mfilename sprintf('Geometric Mean h(x,y)   : %f \n', GM_hxy)]);
LineBreakMedium();
fprintf([mfilename sprintf('Optimal value theta_{1} : %f \n', th1)]);
fprintf([mfilename sprintf('Optimal value theta_{2} : %f \n', th2)]);
fprintf([mfilename sprintf('Optimal value alpha     : %f \n', alpha)]);
fprintf([mfilename sprintf('Optimal value beta      : %f \n', beta)]);
fprintf([mfilename sprintf('Optimal value gamma     : %f \n', gamma)]);
LineBreakLarge();
fprintf([mfilename sprintf(':Degree of GCD : %i \n',t)])
LineBreakLarge();

fxy_n = fxy ./ GM_fxy;
gxy_n = gxy ./ GM_gxy;
hxy_n = hxy ./ GM_hxy;

% %
% Perform low rank approximation method and compute coefficients of the
% polynomials u(x,y) v(x,y) and w(x,y)
[fxy1_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, alpha_lr, ...
    beta_lr, gamma_lr, th1_lr, th2_lr] = GetLowRankApproximation_3Polys(...
    fxy_n, gxy_n, hxy_n, alpha, beta, gamma, th1, th2, m, n, o, t);


% Perform approximate polynomial factorisation to compute the coefficients
% of the GCD d(x,y
[fxy_lra, gxy_lra, hxy_lra, uxy_lra, vxy_lra, wxy_lra, dxy_lra, alpha_lra, beta_lra, gamma_lr, th1_lra, th2_lra] = ...
    APF_3Polys(fxy1_lr, gxy_lr, hxy_lr,...
    uxy_lr, vxy_lr, wxy_lr, ...
    alpha_lr, beta_lr, gamma_lr, th1_lr, th2_lr,...
    m, n, o, t);



% % 
% Remove thetas from u(w,w), v(w,w) and d(w,w) to obtain u(x,y) v(x,y) and
% d(x,y)
fxy_o = fxy_lra;
gxy_o = gxy_lra;
hxy_o = hxy_lra;

dxy_o = dxy_lra;

uxy_o = uxy_lra;
vxy_o = vxy_lra;
wxy_o = wxy_lra;

end

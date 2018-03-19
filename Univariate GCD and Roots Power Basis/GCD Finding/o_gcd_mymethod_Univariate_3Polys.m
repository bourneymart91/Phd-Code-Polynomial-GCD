function [fx_o, gx_o, hx_o, dx_o, ux_o, vx_o, wx_o, alpha_o, theta_o, t , GM_fx, GM_gx, GM_hx] =...
    o_gcd_mymethod_Univariate_3Polys(fx, gx, hx, limits_t, rank_range)
% o_gcd_mymethod(fx,gx,deg_limits)
%
% Given two polynomials f(x) and g(x) return the GCD d(x)
%
% % Inputs.
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% gx : (Vector) Coefficients of polynomial g(x)
%
% hx : (Vector) Coefficients of polynomial h(x)
%
% deg_limits : (Int Int) Specifiy upper and lower bound of the degree of the GCD
% d(x), typically set when o_gcd_mymethod() is called from a root solving
% problem, where deg_limits have been predetermined.
%
% rank_range : [Float Float]
%
%
% % Outputs
% 
% fx_o : (Vector) 
%
% gx_o : (Vector)
%
% dx_o : (Vector)
%
% ux_o : (Vector)
% 
% vx_o : (Vector)
%
% alpha_o : (Float)
%
% theta_o : (Float)




% Preprocess the polynomials f(x) and g(x)
[GM_fx, GM_gx, GM_hx, lambda, mu, rho, theta] = ...
    Preprocess_3Polys(fx, gx, hx);

% Get f(x), g(x) and h(x) normalised by mean
fx_n = fx./ GM_fx;
gx_n = gx./ GM_gx;
hx_n = hx./ GM_hx;

% Get f(\omega), g(\omega) and h(\omega)
fw = GetWithThetas(fx_n, theta);
gw = GetWithThetas(gx_n, theta);
hw = GetWithThetas(hx_n, theta);


% Get degree of polynomials
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

PlotCoefficients({fx, gx, hx, lambda.* fw, mu.*gw, rho.*hw}, ...
    {'$f(x)$','$g(x)$','$h(x)$','$\lambda \tilde{f}(\omega)$',...
    '$\mu\tilde{g}(\omega)$','$\rho\tilde{h}(\omega)$'});

% Get the degree of the GCD
t = GetGCDDegree_3Polys(lambda.*fw, mu.*gw, rho.*hw, limits_t, rank_range);
LineBreakLarge();

% Print the degree of the GCD
fprintf([mfilename ' : ' sprintf('Degree of f(x): % i \n',m)]);
fprintf([mfilename ' : ' sprintf('Degree of g(x): % i \n',n)]);
fprintf([mfilename ' : ' sprintf('Degree of h(x): % i \n',o)]);

fprintf([mfilename ' : ' sprintf('Degree of GCD : % i \n',t)]);



if (t == 0)
    dx_o = 1;
    ux_o = fx;
    vx_o = gx;
    return
end


% %
% %
% Get the low rank approximation and refined values for fx, gx, alpha, and
% theta. Note that vectors of polynomial coefficients are suffixed with 
% 'lr' (lr = low rank) and these are the coefficients after low rank
% approximation is computed.

[fx_lr, gx_lr, hx_lr, ux_lr, vx_lr, wx_lr, alpha_lr, theta_lr] = ...
    LowRankApprox_3Polys(fx_n, gx_n, hx_n, lambda, theta, t);


 [ux_lra, vx_lra, wx_lra, fx_lra, gx_lra, hx_lra, dx_lra, alpha_lra, theta_lra] = ...
     APF_3Polys(ux_lr, vx_lr,wx_lr, fx_lr, gx_lr, hx_lr, alpha_lr, theta_lr, t);

% Get f(x), g(x) and h(x)
fx_o = fx_lra;
gx_o = gx_lra;
hx_o = hx_lra;

% Get u(x), v(x) and w(x)
ux_o = ux_lra;
vx_o = vx_lra;
wx_o = wx_lra;

% Get d(x)
dx_o = dx_lra;

% Get \alpha and \theta
alpha_o = alpha_lra;
theta_o = theta_lra;
end



function [] = PlotCoefficients(arrPolynomials, arrNames)


nPolynomials = length(arrPolynomials);

figure()
hold on
for i = 1:1:nPolynomials
    
    poly = arrPolynomials{i};
    polyName = arrNames{i};
    
    nCoefficients = length(poly);
    x_vec = 1 : 1 : nCoefficients;
    plot(x_vec, log10(abs(poly)), 'DisplayName',polyName);
    
    
end

xlabel('$i$ : Coefficient Index','Interpreter', 'latex', 'FontSize', 20)
ylabel('$\log_{10} \left( \Re \right)$','Interpreter', 'latex', 'FontSize', 20)

l = legend(gca,'show');
set(l,'Interpreter', 'latex')
set(l,'FontSize', 20)

% Figure size and location
myplot = gca;
myval_side = 0.10;
myval_base = 0.08;
set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])
set(gcf, 'Position', [100, 100, 710, 650])

grid on
box on

hold off


end




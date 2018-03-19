function [fxy_o, gxy_o, hxy_o, dxy_o, uxy_o, vxy_o, wxy_o, t, t1, t2, rank_range] = ...
    o_gcd_mymethod_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, limits_t, limits_t1, limits_t2, rank_range)
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
% The polynomials are given matrix form as follows
%
%       1   y   y^2     ...
%        ___________
%   1   |___|___|___|   ...
%   x   |___|___|___|   ...
%   x^2 |___|___|___|   ...
%   ...    .  .  .
%
% % Inputs.
%
%
% fxy : (Matrix) Matrix of the coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Matrix of the coefficients of polynomial g(x,y)
% 
% hxy : (Matrix) Matrix of the coefficients of polynomial h(x,y)
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
% uxy : (Matrix) Calculated coefficients of the polynomial u(x,y) 
%
% vxy : (Matrix) Calculated coefficients of the polynomial v(x,y)
%
% wxy : (Matrix) Calculated coefficients of the polynomial w(x,y)
%
% dxy : Calculated coefficients of d(x,y)

% % Initialise the global variables

% Initialise global settings
global SETTINGS

% % Preprocessing

%[GM_fx, GM_gx, GM_hx, alpha, th1, th2] = Preprocess_Relative_3Polys(fxy,gxy,hxy);
% % Temporary
GM_fx = 1;
GM_gx = 1;
GM_hx = 1;
alpha = 1;
th1 = 1;
th2 = 1;


% Normalise f(x,y) by geometric means
fxy_n = fxy./ GM_fx;
gxy_n = gxy./ GM_gx;
hxy_n = hxy./ GM_hx;

% Get preprocessed forms of f(x,y), g(x,y) and h(x,y)
fww = GetWithThetas(fxy_n, th1, th2);
gww = GetWithThetas(gxy_n, th1, th2);
hww = GetWithThetas(hxy_n, th1, th2);

[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);

limits_k = [0 min([m,n,o])];
limits_k1 = [0 min([m1, n1, o1])];
limits_k2 = [0 min([m2, n2, o2])];



% %
% %
% %
% Given the total degree, compute relative degrees t_{1} and t_{2}

switch SETTINGS.DEGREE_METHOD
    case 'Total'
        % Must compute t1 and t2 so do so with relative method
        [t, rank_range] = GetGCDDegree_Total_Bivariate_3Polys(fww, alpha.*gww, hww, ...
            m, n, o, limits_t, rank_range);
        t1 = 1000;
        t2 = 1000;
        
        %[t1,t2] = GetGCDDegree_Relative(fww,alpha*gww,m,n,t);
        
    case 'Relative'
        
        t = 1000;
        [t1, t2] = GetGCDDegree_Relative_Bivariate_3Polys(fww, alpha.*gww, ...
            hww, limits_k1, limits_k2, limits_t1, limits_t2, rank_range);
        
    case 'Both'
        
        [t] = GetGCDDegree_Total_Bivariate_3Polys(fww, alpha.*gww, hww,...
            m, n, o, limits_t, rank_range);
        
        [t1, t2] = GetGCDDegree_Relative_Given_t_Bivariate_3Polys(fww, alpha*gww, hww, m , n, o, t, limits_t1, limits_t2, rank_range);
        
    otherwise
        error('%s is not a valid DEGREE_METHOD', SETTINGS.DEGREE_METHOD)
end



LineBreakLarge();
fprintf([mfilename ' : ' sprintf('The calculated relative degree is : t_{1} = %i, t_{2} = %i \n',t1,t2)])
fprintf([mfilename ' : ' sprintf('The calculated total degree is :  %i \n',t)])
LineBreakLarge();



% %
% %
% %
% Get optimal column for removal from S_{t_{1},t_{2}}
idx_col = GetOptimalColumn_3Polys(fww, alpha.*gww, hww, m, n, o, t, t1, t2);


% % 
% %
% %
% Perform low rank approximation method
[fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, alpha_lr, th1_lr, th2_lr] = ...
    LowRankApprox_Bivariate_3Polys(fxy_n, gxy_n, hxy_n, alpha, th1, th2, m, n, o, t, t1, t2, idx_col);



[fxy_lra, gxy_lra, hxy_lra, uxy_lra, vxy_lra, wxy_lra, dxy_lra, alpha_lra, th1_lra, th2_lra] = ...
    APF_Bivariate_3Polys(fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, m, n, o, t, alpha_lr, th1_lr, th2_lr);



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

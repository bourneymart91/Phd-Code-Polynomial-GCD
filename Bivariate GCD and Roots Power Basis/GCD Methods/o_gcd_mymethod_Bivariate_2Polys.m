function [fxy_o, gxy_o, dxy_o, uxy_o, vxy_o, t, t1, t2] = o_gcd_mymethod_Bivariate_2Polys(fxy, gxy, m, n, limits_t, limits_t1, limits_t2, rank_range)
% o_gcd_mymethod_Bivariate_2Polys(fxy, gxy, m, n, myLimits_t, myLimits_t1, myLimits_t2, limits_t, limits_t1, limits_t2)
%
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
% fxy : (Matrix) of coefficients of polynomial f(x,y)
%
% gxy : (Matrix) of coefficients of polynomial g(x,y)
%
% m : (Int) Total degree of polynomial f(x,y)
%
% n : (Int) Total degree of polynomial f(x,y) and g(x,y)
%
% limits_t : (Int Int)
%
% limits_t1 : (Int Int)
%
% limits_t2 : (Int Int)
% 
% % Outputs
%
% uxy : (Matrix) Coefficients of the polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of the polynomial v(x,y)
%
% dxy : (Matrix) Coefficients of the polynomial d(x,y)




% % Initialise the global variables

% Initialise global settings
global SETTINGS

% % Preprocessing

[GM_fxy, GM_gxy, alpha, th1, th2] = Preprocess_Relative(fxy, gxy);

%[lambda_b, mu_b, alpha_b, th1_b th2_b] = Preprocess_Total(fxy, gxy, m, n);
%
%[lambda_c, mu_c, alpha_c, th1_c th2_c] = Preprocess_Both(fxy, gxy, m, n);

% Normalise f(x,y) by geometric means
fxy_n = fxy ./ GM_fxy;
gxy_n = gxy ./ GM_gxy;

% Get f(w,w) from f(x,y) and g(w,w) from g(x,y)
fww = GetWithThetas(fxy_n, th1, th2);
gww = GetWithThetas(gxy_n, th1, th2);


switch SETTINGS.DEGREE_METHOD
    case 'Total'
        
        % Compute the total degree of the GCD
        [t] = GetGCDDegree_Total_Bivariate_2Polys(fww, alpha.*gww, m, n, limits_t, rank_range);
        [t] = GetGCDDegree_Total_Bivariate_2Polys_Fast(fww, alpha.*gww, m, n, limits_t, rank_range);
        
        % set relative degrees of GCD to arbitrary values
        t1 = 1000;
        t2 = 1000;
        
    case 'Relative'
        
        % Set total degree to an arbitrary value
        t = 1000;
        
        % Compute the relative degree of the GCD
        [m1, m2] = GetDegree_Bivariate(fxy);
        [n1, n2] = GetDegree_Bivariate(gxy);
        
        
        limits_t1 = [0 min(m1, n1)];
        limits_t2 = [0 min(m2, n2)];
        %
        %tic;
        [t1 ,t2] = GetGCDDegree_Relative_Bivariate_2Polys(fww, alpha.*gww, limits_t1, limits_t2);
        %toc;
%         % NEW 15/02/2017
        %tic;
        [t1, t2] = GetGCDDegree_Relative_Bivariate_2Polys_Fast(fww, alpha.*gww, limits_t1, limits_t2);
        %toc;
      
        
    case 'Both'
        
        % Compute total degree
        [t] = GetGCDDegree_Total_Bivariate_2Polys(fww, alpha.*gww, m, n, limits_t, rank_range);
        
        % Compute relative degree given total degree
        
        % Set degree Limits
        [m1, m2] = GetDegree_Bivariate(fxy);
        [n1, n2] = GetDegree_Bivariate(gxy);
        
        limits_t1 = [0 min(m1,n1)];
        limits_t2 = [0 min(m2,n2)];
        
        % Compute relative degree given total degree
        [t1, t2] = GetGCDDegree_Relative_Given_t_Bivariate_2Polys(fww, alpha*gww, m, n, t, limits_t1, limits_t2);
        
        %[t1_test,t2_test] = GetGCDDegree_Relative_Bivariate_2Polys(fww, alpha*gww, limits_t1, limits_t2);
        
    otherwise
        error('err')
end



LineBreakMedium();
fprintf([mfilename ' : ' sprintf('The calculated relative degree is : t_{1} = %i, t_{2} = %i \n',t1, t2)])
LineBreakMedium();



% %
% %
% %
% Get optimal column for removal from Sylvester matrix S_{t} or S_{t1,t2}
% or S_{t,t1,t2}
idx_col = GetOptimalColumn_2Polys(fww, alpha.*gww, m, n, t, t1, t2);


% % 
% %
% %
% Perform low rank approximation method
[fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = LowRankApprox_Bivariate_2Polys(fxy_n, gxy_n, alpha, th1, th2, m, n, t, t1, t2, idx_col);



[fxy_lra, gxy_lra, uxy_lra, vxy_lra, dxy_lra, alpha_lra, th1_lra, th2_lra] = ...
    APF_Bivariate_2Polys(fxy_lr, gxy_lr, uxy_lr, vxy_lr, m, n, t, alpha_lr, th1_lr, th2_lr);



% % 
% Remove thetas from u(w,w), v(w,w) and d(w,w) to obtain u(x,y) v(x,y) and
% d(x,y)
fxy_o = fxy_lra;
gxy_o = gxy_lra;
dxy_o = dxy_lra;
uxy_o = uxy_lra;
vxy_o = vxy_lra;


end

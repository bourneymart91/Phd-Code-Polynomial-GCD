function [root_mult_array] = o_roots_Yun(fx)
% Given the polynomial f(x) compute its roots by Square free factorization.
% This algorithm is referred to as Yuns Algorithm in
% https://en.wikipedia.org/wiki/Square-free_polynomial
%
%
% % Inputs
%
% fx : Coefficients of polynomial f(x) as a column vector where first
% coefficient has lowest power [a_{0} ; ... ; a_{m}]^{T}
%
% % Outputs
%
% root_mult_array : A matrix of roots and multiplicities, where the first
% column contains all of the roots of f(x) and the second contains the
% corresponding multiplicities.

global SETTINGS
SETTINGS.PLOT_GRAPHS = true;

% Set Iteration number
ite = 1;

f_dash = Differentiate(fx);

% Get degree of f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(f_dash);

limits_k = [1,min(m,n)];

% Perform GCD computation.

[fx_n,gx_n,uxy, C{ite}, arr_wx{ite}, alpha, theta, t , lambda,mu] ...
    = o_gcd_mymethod_Univariate_2Polys(fx, Differentiate(fx), limits_k);

arr_dxy{ite} = arr_wx{ite} - Differentiate(C{ite});
LineBreakMedium();


while (GetDegree(C{ite}) > 0 )
    
    % Get the degree of polynomial f(x)
    m = GetDegree(C{ite});
    % Get the degree of polynomial g(x)
    n = GetDegree(arr_dxy{ite});
    
    limits_t = [1,min(m,n)];
    
    % Get GCD
    [fx,gx, arr_hx{ite+1},C{ite+1},arr_wx{ite+1},~,~,~,~,~] =...
        o_gcd_mymethod_Univariate_2Polys(C{ite}, arr_dxy{ite}, limits_t);
    
   
    
    arr_dxy{ite+1} = arr_wx{ite} - Differentiate(C{ite});
    
    ite = ite+1;
    
    LineBreakMedium();
    
end

SETTINGS.PLOT_GRAPHS = false;

root_mult_array = [];

for i = 1:1:length(arr_hx)

    try
        
    factor = arr_hx{i};
    % Divide by x coefficient
    factor = factor./factor(2);
    % Get root
    root = -1.*factor(1)
    
    catch
        
    end
end


end
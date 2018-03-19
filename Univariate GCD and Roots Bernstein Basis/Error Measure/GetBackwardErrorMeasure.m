
function [my_error] = GetBackwardErrorMeasure(root_mult_arr, fx_exact)
% Get the distance between the polynomial f_{comp} and f_{exact}
%
% Inputs.
%
% root_mult_arr : Matrix whose rows contain a computed root of f(x) and
%                 its corresponding multiplicitiy.
%
% fx_exact : Coefficients of exact Polynomial f(x)
%
% % Outputs
%
% Error : Measure of error between the exact polynomial f(x) and the
% polynomial f(x) computed by the root finding method.

% Get coefficients of computed polynomial f(x)
fx_comp = GetWithoutBinomials( BuildPolyFromRoots( root_mult_arr));

% Get distance between f_{comp}(x) and f_{exact}(x)

fx_comp = fx_comp./fx_comp(1);
fx_exact = fx_exact ./ fx_exact(1);

my_error  = norm(fx_comp - fx_exact) ./ norm(fx_exact);


end
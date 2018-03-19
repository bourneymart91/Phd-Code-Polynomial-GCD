function [myError] = GetPolynomialError(fx_exact, fx_calc)
% Used in Computing the distance between the GCD and the computed GCD in
% the GCD finding method.
%
% Get distance between f(x) and the calulated f(x)
%
% % Inputs
%
% name : (String) Name of the polynomial eg 'f(x)'
%
% fx_calc : (Vector) The coefficients of the polynomial f(x) as calculated
%
% fx_exact : (Vector) The coefficients of the exact polynomial f(x)


% Check both vectors are column vectors of the same size
[nRows1, nCols1] = size(fx_exact);
[nRows2, nCols2] = size(fx_calc);

if (nRows1 ~= nRows2) || (nCols1 ~= nCols2)
   myError = 10000000000;
   display('Error in Comparing Polynomials : GetPolynomialError()')
   return
end



% Get the angle between the two vectors
% angle = dot(f_calc,f_exact) ./ (norm(f_calc) * norm(f_exact));
% angle_error = 1 - angle;
% fprintf('\tCalculated angle error : %8.2e \n', angle_error)

% Get unit vectors
fx_calc  = NormaliseVector(fx_calc);
fx_exact = NormaliseVector(fx_exact);


% Make leading coefficient positive for both vectors
fx_calc = MakePositive(fx_calc);
fx_exact = MakePositive(fx_exact);

try
    myError = norm(abs(fx_exact - fx_calc)) ./ norm(fx_exact);
catch
    myError = 1000;
end





end



function fx = MakePositive(fx)

if fx(1) < 0
    fx = fx.*-1;
end

end
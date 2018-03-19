function [] = o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_variant,  rank_revealing_metric)
% o_gcd_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
%   low_rank_approx_method, apf_method, sylvester_matrix_variant, ...
%   rank_revealing_metric)
%
% Obtain the Greatest Common Divisor (GCD) d(x) of two polynomials f(x) and
% g(x) as defined in the example file.
%
%
% % Inputs.
%
% ex_num : (String) Example Number
%
% emin: (Float) Signal to noise ratio (minimum)
%
% emax: (Float) Signal to noise ratio (maximum)
%
% mean_method : (String) Method for taking mean of entires in S_{k}
%   * 'Geometric Mean Matlab Method'
%   * 'Geometric Mean My Method'
%
% bool_alpha_theta : (Bool) true or false if preprocessing is performed
%   * true
%   * false
%
% low_rank_approx_method : (String)
%   * 'Standard STLN'
%   * 'Standard SNTLN'
%   * 'Root Specific SNTLN'
%   * 'None'
%
% apf_method : (String)
%   * 'Standard APF NonLinear'
%   * 'Standard APF Linear'
%   * 'None'
%
% sylvester_matrix_variant : (String)
%   * 'T'
%   * 'DT'
%   * 'DTQ'
%   * 'TQ'
%   * 'DTQ Denominator Removed'
%   * 'DTQ Rearranged'
%
% nEquations : (String)
%   * '2' : 2 Equations define the 2 x 3 Sylvester matrix
%   * '3' : 3 Equations define the 3 x 3 Sylvester matrix
%
% rank_revealing_metric : (String)
%   * 'Minimum Singular Values' :
%   * 'Max/Min Singular Values' :
%   * 'R1 Row Norms' :
%   * 'R1 Row Diagonals' :
%   * 'Residuals' :
%
% % Example
% >> o_gcd_Univariate_2Polys('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'None', 'None', 'DTQ', 'Minimum Singular Values')
% >> o_gcd_Univariate_2Polys('1',1e-10,1e-12,'Geometric Mean Matlab Method', true, 'Standard STLN', 'Standard APF Nonlinear', 'DTQ', 'Minimum Singular Values')

global SETTINGS

% Add that folder plus all subfolders to the path.
addpath(genpath(pwd));

% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end


% % Set global variables.
SetGlobalVariables_GCD_2Polys(...
    ex_num,...
    emin,...
    emax,...
    mean_method,...
    bool_alpha_theta,...
    low_rank_approx_method,...
    apf_method,...
    sylvester_matrix_variant, ... 
    rank_revealing_metric);

% Print the parameters.
LineBreakLarge()
fprintf('INPUT VARIABLES\n')
fprintf('\t EXAMPLE NUMBER : %s \n',ex_num)
fprintf('\t EMIN : %s \n' , num2str(SETTINGS.EMIN))
fprintf('\t EMAX : %s \n' , num2str(SETTINGS.EMAX))
fprintf('\t MEAN METHOD : %s \n', SETTINGS.MEAN_METHOD)
fprintf('\t ALPHA_THETA : %s \n', num2str(SETTINGS.BOOL_ALPHA_THETA))
fprintf('\t RANK REVEALING METRIC : %s \n', SETTINGS.RANK_REVEALING_METRIC);
fprintf('\t LOW RANK APPROX METHOD : %s \n', SETTINGS.LOW_RANK_APPROXIMATION_METHOD);
fprintf('\t APF METHOD : %s \n ', SETTINGS.APF_METHOD);
fprintf('\t LOG: %s \n', num2str(SETTINGS.BOOL_LOG));
fprintf('\t SYLVESTER MATRIX VARIANT : %s \n', SETTINGS.SYLVESTER_MATRIX_VARIANT);
LineBreakLarge()

% o - gcd - Calculate GCD of two Arbitrary polynomials
% Given two sets of polynomial roots, form polynomials f and g, expressed
% in the Bernstein Basis. Add noise, and calculate the GCD of the two
% polynomails



% Get roots from example file
[fx_exact, gx_exact, dx_exact, ux_exact, vx_exact] = Examples_GCD(ex_num);






% Get degree of f(x), g(x) and d(x)
m = GetDegree(fx_exact);
n = GetDegree(gx_exact);
t_exact = GetDegree(dx_exact);

% Add componentwise noise to coefficients of polynomials in 'Standard Bernstein Basis'.
fx_noisy = AddVariableNoiseToPoly(fx_exact, emin, emax);
gx_noisy = AddVariableNoiseToPoly(gx_exact, emin, emax);

% Set upper and lower limit of the degree of the GCD
lowerLimit_t = 1;
upperLimit_t = min(m, n);
limits_t = [lowerLimit_t upperLimit_t];

% Set range of rank revealing metric
rank_range = [0 0];




% Obtain the coefficients of the GCD d and quotient polynomials u(x) and v(x).
[fx_calc, gx_calc, dx_calc, ux_calc, vx_calc, alpha_calc, theta_calc, t_calc, rank_range] = ...
    o_gcd_2Polys_mymethod(fx_noisy, gx_noisy, limits_t, rank_range);


%
% BEZOUT METHOD
%

if SETTINGS.BOOL_BEZOUTIAN

    bool_preproc = true;
    [dx_calc_bez_qr, dx_calc_bez_lu] = o_mod(fx_noisy, gx_noisy, bool_preproc);
    error.dx_bez_qr = GetPolynomialError(dx_exact, (dx_calc_bez_qr));
    error.dx_bez_lu = GetPolynomialError(dx_exact, (dx_calc_bez_lu));

    display(error.dx_bez_qr)
    display(error.dx_bez_lu)


end


%[d, u, v] = o_GCD_matlab(fx_noisy, gx_noisy);


% Check coefficients of calculated polynomials are similar to those of the
% exact polynomials.

LineBreakMedium();
try
    
    % Get error between \hat{f}(x) and inexact polynomial f(x)
    error.fx = GetPolynomialError(fx_exact, fx_noisy);
    error.fx_calc = GetPolynomialError(fx_exact, fx_calc);
    
    % Get error between \hat{g}(x) and inexact polynomial g(x)
    error.gx = GetPolynomialError(gx_exact, gx_noisy);
    error.gx_calc = GetPolynomialError(gx_exact, gx_calc);
    
    %
    error.ux = GetPolynomialError(ux_exact, ux_calc);
    error.vx = GetPolynomialError(vx_exact, vx_calc);
    error.dx = GetPolynomialError(dx_exact, dx_calc);
    
    %
    error.uw = GetPolynomialError(GetWithThetas(ux_exact, theta_calc), GetWithThetas(ux_calc, theta_calc));
    error.vw = GetPolynomialError(alpha_calc * GetWithThetas(vx_exact, theta_calc), alpha_calc * GetWithThetas(vx_calc, theta_calc));
    error.dw = GetPolynomialError(GetWithThetas(dx_exact, theta_calc), GetWithThetas(dx_calc, theta_calc));
    
    LineBreakMedium()
    fprintf('Error u(x) : %e \n', error.ux);
    fprintf('Error v(x) : %e \n', error.vx);
    fprintf('Error d(x) : %e \n', error.dx);
    fprintf('Average Error : %e \n', mean([error.ux, error.vx, error.dx]))
    LineBreakMedium()
    fprintf('Error u(omega) : %e \n', error.uw);
    fprintf('Error v(omega) : %e \n', error.vw);
    fprintf('Error d(omega) : %e \n', error.dw);
    fprintf('Average Error : %e \n', mean([error.uw, error.vw, error.dw]))
    LineBreakMedium()
    
catch
    
    fprintf('Error : Can not compare computed GCD with Exact GCD. Check degree of GCD is computed correctly \n')
    error.dx = 1000;
    error.ux = 1000;
    error.vx = 1000;
    
end

% Print results to results file
PrintGCDToFile(m, n, t_exact, t_calc, error);

LineBreakMedium();

end






function [] = CompareCoefficients(fx_exact, fx_noisy, fx_calc, str1, str2, str3) 


fx_exact_unit = fx_exact ./ norm(fx_exact);
fx_noisy_unit = fx_noisy ./ norm(fx_noisy);
fx_calc_unit = fx_calc ./ norm(fx_calc);


PlotCoefficients(...
    {...
    fx_exact_unit, ...
    fx_noisy_unit, ...
    fx_calc_unit...
    abs(fx_exact_unit - fx_calc_unit)./fx_exact_unit...
    abs(fx_exact_unit - fx_calc_unit)./fx_exact_unit...
    },...
    {...
    str1,...
    str2,...
    str3,...
    sprintf('distance %s and %s',str1,str2)...
    sprintf('distance %s and %s',str1,str3)...
    }...
    ,...
    {...
      '-o', ...
      '-s', ...
      '-', ...
      '-', ...
      '-'...
    }...
);


end




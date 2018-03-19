function [] = o_gcd_Univariate_3Polys(ex_num, ex_num_variant, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_variant, nEquations, rank_revealing_metric)
% o_gcd_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
%   low_rank_approx_method, apf_method)
%
% Obtain the Greatest Common Divisor (GCD) d(x) of two polynomials f(x) and
% g(x) as defined in the example file.
%
%
% % Inputs.
%
% ex : (String) Example Number
%
% ex_num_variant : (String) 'a', 'b' or 'c' determines the ordering of the
%   polynomials f(x), g(x) and h(x) in the subresultant matrices.
%
% emin: (Float) Signal to noise ratio (minimum)
%
% emax: (Float) Signal to noise ratio (maximum)
%
% mean_method : Method for taking mean of non-zero entires in the 
%   subresutlant matrices S_{k}.
%   'Geometric Mean Matlab Method'
%   'Geometric Mean My Method'
%   'None'
%
% bool_alpha_theta : (Bool) Boolean value determines whether the
%   subresultant matrices of the subresultant matrix sequence are
%   preprocessed. Note the variable name refers to alpha and theta, but
%   optimal values of \lamda_{k}, \mu_{k}, \rho_{k} and \theta_{k} are
%   computed. 
%           true : Subresultant matrices are to be preprocessed
%           false : Subresultant matrices are not preprocessed
%
% low_rank_approx_method : (String)
%   'Standard STLN'
%   'Standard SNTLN'
%   'Root Specific SNTLN'
%   'None'
%
% apf_method : (String)
%   'Standard APF NonLinear'
%   'Standard APF Linear'
%   'None'
%
% sylvester_matrix_variant : (String)
%   'T'
%   'DT'
%   'DTQ'
%   'TQ'
%   'DTQ Denominator Removed'
%   'DTQ Rearranged'
%
% nEquations : (String) The number of equations used to define the
%   subresutlant matrices.
%   '2' : Syvlester matrix defined by 2 equations, and has 2 x 3 structure
%   '3' : Sylvester matrix defined by 3 equations, and has 3 x 3 structure
%
% rank_revealing_metric : (String)
%   'Minimum Singular Values'
%
%
%
%
% % Example
%
%
% >> o_gcd_Univariate_3Polys('1', 'a', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'None', 'None', 'DTQ', '2', 'Minimum Singular Values')
%
% >> o_gcd_Univariate_3Polys('1', 'a', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'Standard STLN', 'Standard APF Nonlinear', 'DTQ', '2', 'Minimum Singular Values')
%


% if not 11 arguments then error
if (nargin ~= 11)
    error('Not enough input arguments')
end



% Add subfolders
restoredefaultpath
addpath(genpath(pwd))


% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end

% % Set global variables.
SetGlobalVariables_GCD_3Polys(...
    ex_num,...
    ex_num_variant,...
    emin,...
    emax,...
    mean_method,...
    bool_alpha_theta,...
    low_rank_approx_method,...
    apf_method,...
    sylvester_matrix_variant, ...
    nEquations,...
    rank_revealing_metric);




% o - gcd - Calculate GCD of two Arbitrary polynomials
% Given two sets of polynomial roots, form polynomials f and g, expressed
% in the Bernstein Basis. Add noise, and calculate the GCD of the two
% polynomails



% Get the polynomials f(x), g(x), h(x), the GCD d(x), and the cofactor 
% polynomials u(x), v(x) and w(x) from the example file.
[fx_exact, gx_exact, hx_exact, dx_exact, ux_exact, vx_exact, wx_exact] = ...
    Examples_GCD_3Polys(ex_num, ex_num_variant);



% Add componentwise noise to coefficients of the polynomials f(x), g(x) and
% h(x).
fx = AddVariableNoiseToPoly(fx_exact, emin, emax);
gx = AddVariableNoiseToPoly(gx_exact, emin, emax);
hx = AddVariableNoiseToPoly(hx_exact, emin, emax);


%PlotCoefficients({fx, gx, hx}, {'fx','gx','hx'});

% Get the degree of the polynomials f(x), g(x) and h(x) and the degree of
% the exact GCD d(x).
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);
t_exact = GetDegree(dx_exact);


LineBreakLarge()
fprintf('\t m : %i \n', m)
fprintf('\t n : %i \n', n)
fprintf('\t o : %i \n', o)
fprintf('\t t : %i \n', t_exact)
LineBreakLarge()





% Set the upper and lower limit of the degree of the GCD. Since this is a
% general GCD problem, no prior limits are known so these values are set to
% 1 and min(m,n,o).
lower_limit_t = 1;
upper_limit_t = min([m n o]);
limits_t = [lower_limit_t, upper_limit_t];



% Set the rank_range. Again, since this is a stand alone GCD problem limits
% of the rank_revealing_metric are not known, so set to [0,0]
rank_range = [0 0];






% Obtain the coefficients of the GCD d and quotient polynomials u and v.
[fx_calc, gx_calc, wx_calc, dx_calc, ux_calc, vx_calc, wx_calc, ...
    lambda_calc, mu_calc, rho_calc, theta_calc] = ...
    o_gcd_3Polys_mymethod(fx, gx, hx, limits_t, rank_range);





% Check coefficients of calculated polynomials are similar to those of the
% exact polynomials.



my_error.dx = GetPolynomialError(dx_exact, dx_calc);
my_error.ux = GetPolynomialError(ux_exact, ux_calc);
my_error.vx = GetPolynomialError(vx_exact, vx_calc);
my_error.wx = GetPolynomialError(wx_exact, wx_calc);


LineBreakMedium();
fprintf('Error u(x) : %e \n', my_error.ux);
fprintf('Error v(x) : %e \n', my_error.vx);
fprintf('Error w(x) : %e \n', my_error.wx);
fprintf('Error d(x) : %e \n', my_error.dx);
fprintf('Average Error : %e \n', mean([my_error.ux, my_error.vx, my_error.dx]))
LineBreakMedium();


% Print results to results file
PrintResultToFile(GetDegree(fx), GetDegree(gx), GetDegree(hx), GetDegree(dx_calc), my_error)



end




function [] = PrintResultToFile(m, n, o, t, my_error)
% Print Results to file
%
% % Inputs
%
% m : (Int) The degree of the polynomial f(x)
%
% n : (Int) The degree of the polynomial g(x)
%
% o : (Int) The degree of the polynomial h(x)
%
% t : (Int) The computed degree of the GCD d(x)
%
% error : array of errors e
%   error.dx : (Float)
%   error.ux : (Float)
%   error.vx : (Float)
%   error.wx : (Float)


global SETTINGS

% Specify file name
fullFileName = sprintf('Results/Results_o_gcd_3Polys.dat');

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName, 'a');
    WriteNewLine()
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    WriteNewLine()
    fclose(fileID);
    
end



    function WriteNewLine()
        
        fprintf(fileID,'%s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s \n',...
            datetime('now'),...
            SETTINGS.EX_NUM,...
            int2str(m),...
            int2str(n),...
            int2str(o),...
            int2str(t),...
            num2str(my_error.ux),...
            num2str(my_error.vx),...
            num2str(my_error.wx),...
            num2str(my_error.dx),...
            SETTINGS.MEAN_METHOD,...
            num2str(SETTINGS.BOOL_ALPHA_THETA),...
            SETTINGS.EMIN,...
            SETTINGS.EMAX,...
            SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
            num2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
            SETTINGS.APF_METHOD,...
            num2str(SETTINGS.APF_REQ_ITE),...
            SETTINGS.BOOL_LOG,...
            SETTINGS.SYLVESTER_MATRIX_VARIANT,...
            SETTINGS.GCD_COEFFICIENT_METHOD...
            );
        % 21 inputs
        
    end

    function WriteHeader()
        strHeaders = ['DATE, EX_NUM, m, n, o, t, ERROR_UX, ERROR_VX,' ...
            'ERROW_WX, ERROR_DX, MEAN_METHOD, BOOL_ALPHA_THETA, EMIN, '...
            'EMAX, LOW_RANK_APPROX_METHOD, LOW_RANK_ITE, APF_METHOD, '...
            'APF_ITE, BOOL_LOG, SYLVESTER_MATRIX_VARIANT, GCD_METHOD\n'...
        ];
    
        fprintf(fileID, strHeaders);
    end


end


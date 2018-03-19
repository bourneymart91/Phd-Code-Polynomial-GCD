function [] = o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method,...
    bool_alpha_theta, low_rank_approx_method, apf_method, degree_method, ...
    rank_revealing_metric)
% o_gcd_bivar_2Polys(ex_num,emin,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method,degree_method)
%
% Calculate the GCD d(x,y) of two polynomials f(x,y) and g(x,y) taken from
% the example file.
%
% % Inputs.
%
% ex_num  : (String) Example Number (String)
%
% emin : (Float) Minimum Noise level
%
% emax : (Float) Maximum signal to noise ratio
%
% mean_method : (String)
%   * 'Geometric Mean Matlab Method'
%   * 'None'
%
% bool_alpha_theta (Boolean)
%   * true : Include Preprocessing
%   * false : Exclude Preprocessing
%
% low_rank_approx_method (String)
%   * 'Standard SNTLN'
%   * 'None'
%
% apf_method (String)
%   * 'None'
%   * 'Standard APF Linear'
%   * 'Standard APF Nonlinear'
%
%
% degree_method (String)
%   * 'Relative' : Define polynomials in terms of degree with respect to
%                  x and y, so matrices of coefficients are rectangular.
%   * 'Total' : Define polynomials in terms of total degree, so matrices
%               of coefficients are square, where lower right triangle
%               are all zero.
%   * 'Both' :  Combination of both of the above, and typically gives best
%               results.
%
%
% rank_revealing_metric (String)
%
%
% % Example
% >> o_gcd_Bivariate_2Polys('1',1e-10,1e-12, 'Geometric Mean Matlab Method', true, 'None', 'None', 'Both')
% >> o_gcd_Bivariate_2Polys('1',1e-10,1e-12, 'Geometric Mean Matlab Method', true, 'Standard STLN', 'Standard Nonlinear APF', 'Relative')


% Set the Global Variables
global SETTINGS

% Add subfolders
restoredefaultpath

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 

% Add that folder plus all subfolders to the path.
addpath(genpath(folder));


problem_type = 'GCD';

% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end

% Set global variables
SetGlobalVariables_GCD(problem_type, ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, degree_method,...
    rank_revealing_metric)

% Print Parameters to screen
fprintf('INPUTS. \n')
fprintf('EXAMPLE NUMBER %s \n',ex_num)
fprintf('EMIN : %s \n',emin)
fprintf('EMAX : %s \n',emax)
fprintf('MEAN METHOD : %s \n', mean_method)
fprintf('PREPROCESSING : %s \n',int2str(bool_alpha_theta))
fprintf('RANK_REVEALING_METRIC: %s \n', rank_revealing_metric)
fprintf('LOW RANK METHOD : %s \n',low_rank_approx_method)
fprintf('APF METHOD : %s \n', apf_method)
fprintf('DEGREE METHOD : %s \n', degree_method)

% Get example polynomials
[fxy_exact, gxy_exact,...
    uxy_exact, vxy_exact,...
    dxy_exact,...
    m, m1, m2,...
    n, n1, n2,...
    t_exact, t1_exact, t2_exact] = Examples_GCD_Bivariate_2Polys(ex_num);


DisplayDegreeStructure();

% %
% %
% Add Noise

% Add noise to the coefficients of f and g
[fxy, ~] = AddVariableNoiseToPoly(fxy_exact, emin, emax);
[gxy, ~] = AddVariableNoiseToPoly(gxy_exact, emin, emax);

% %
% % Get the GCD by zengs method
%[u,v,w] = o_gcd_zeng(fxy,gxy);



switch SETTINGS.DEGREE_METHOD
    case 'Total'
        
        fxy_matrix_padd = zeros(m + 1, m + 1);
        gxy_matrix_padd = zeros(n + 1, n + 1);
        
        [nRows_f,nCols_f] = size(fxy);
        fxy_matrix_padd(1 : nRows_f, 1 : nCols_f) = fxy;
        
        [nRows_g,nCols_g] = size(gxy);
        gxy_matrix_padd(1 : nRows_g, 1 : nCols_g) = gxy;
        
        fxy = fxy_matrix_padd;
        gxy = gxy_matrix_padd;
        
    case 'Relative'
        
    case 'Both'
        
    otherwise
        error('err')
end



% Set Limits for the values of k 
% Difference between myLimits and 'limits'. My Limits can be redefined,
% limits should always be computed by number of distinct roots rule.
% Since this is a GCD problem, set 'limits' to default 0,...,min(m,n)
limits_t = [0 min(m, n)];
limits_t1 = [0 min(m1, n1)];
limits_t2 = [0 min(m2, n2)];


rank_range_low = 0;
rank_range_high = 0;
rank_range = [rank_range_low rank_range_high];


% Get the GCD by my method
[fxy_calc, gxy_calc, dxy_calc, uxy_calc, vxy_calc, t, t1, t2] = ...
    o_gcd_mymethod_Bivariate_2Polys(fxy, gxy, m, n, limits_t, limits_t1, limits_t2, rank_range );




% Compare exact d(x,y) and calculated d(x,y)
% PrintCoefficients(dxy_calc,dxy_exact,'d(x,y)');

switch SETTINGS.DEGREE_METHOD
    case 'Relative'
        
        LineBreakLarge()
        my_error.dxy = GetDistanceBetweenPolynomials(dxy_exact, dxy_calc, 'd(x,y)');
        my_error.uxy = GetDistanceBetweenPolynomials(uxy_exact, uxy_calc, 'u(x,y)');
        my_error.vxy = GetDistanceBetweenPolynomials(vxy_exact, vxy_calc, 'v(x,y)');
        LineBreakLarge()
        
    case 'Total'
        
        % Get d(x,y) in a matrix in terms of total degree.
        dxy_exact_total = zeros(t_exact+1,t_exact+1);
        dxy_exact_total(1:t1_exact+1,1:t2_exact+1) = dxy_exact;
        
        uxy_exact_total = zeros(m-t_exact+1, m-t_exact+1);
        uxy_exact_total(1:m1-t1_exact+1,1:m2-t2_exact+1) = uxy_exact;
        
        vxy_exact_total = zeros(n-t_exact+1, n-t_exact+1);
        vxy_exact_total(1:n1-t1_exact+1,1:n2-t2_exact+1) = vxy_exact;
        
        LineBreakLarge()
        my_error.dxy = GetDistanceBetweenPolynomials(dxy_exact_total, dxy_calc, 'd(x,y)');
        my_error.uxy = GetDistanceBetweenPolynomials(uxy_exact_total, uxy_calc, 'u(x,y)');
        my_error.vxy = GetDistanceBetweenPolynomials(vxy_exact_total, vxy_calc, 'v(x,y)');
        LineBreakLarge()
        
    case 'Both'
        
        LineBreakLarge()
        my_error.dxy = GetDistanceBetweenPolynomials(dxy_exact, dxy_calc, 'd(x,y)');
        my_error.uxy = GetDistanceBetweenPolynomials(uxy_exact, uxy_calc, 'u(x,y)');
        my_error.vxy = GetDistanceBetweenPolynomials(vxy_exact, vxy_calc, 'v(x,y)');
        LineBreakLarge()
        
end



PrintToFile(m, n, t, t1, t2, my_error)



end





function []= PrintToFile(m, n, t, t1, t2, error)

global SETTINGS

fullFileName = sprintf('Results/Results_o_gcd.txt');

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    WriteNewLine()
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    WriteNewLine()
    fclose(fileID);
    
end



    function WriteNewLine()
       fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
        datetime('now'),...
        SETTINGS.EX_NUM,...
        int2str(m),...
        int2str(n),...
        int2str(t),...
        int2str(t1),...
        int2str(t2),...
        num2str(error.uxy),...
        num2str(error.vxy),...
        num2str(error.dxy),...
        SETTINGS.MEAN_METHOD,...
        num2str(SETTINGS.BOOL_ALPHA_THETA),...
        SETTINGS.EMIN,...
        SETTINGS.EMAX,...
        SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
        int2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
        SETTINGS.APF_METHOD,...
        SETTINGS.DEGREE_METHOD, ...
        SETTINGS.RANK_REVEALING_METRIC); 
    end


    function WriteHeader()
        fprintf(fileID,'Date, EX_NUM, m, n, t, t1, t2, error_uxy, error_vxy, error_dxy, mean_method, bool_alpha_theta, emax, emin, Low Rank Approx Method, Ite, APF Method, Degree Method, Rank Revealing Metric\n');
    end








end
function [] = o_gcd_Bivariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, degree_method)
% o_gcd(ex_num,emin,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method,degree_method)
%
% Calculate the GCD d(x,y) of two polynomials f(x,y) and g(x,y) taken from
% the example file.
%
% % Inputs.
%
% ex_num : (String) Example Number (String)
%
% emin : (Float) Minimum Noise level
%
% emax : (Float) Maximum signal to noise ratio
%
% mean_method : (String)
%       'Geometric Mean Matlab Method'
%       'None'
%
% bool_alpha_theta (Boolean)
%       true - Include Preprocessing
%       false - Exclude Preprocessing
%
% low_rank_approx_method (String)
%       'Standard SNTLN'
%       'None'
%
% apf_method (String)
%       'None'
%       'Standard APF Linear'
%       'Standard APF Nonlinear'
%
%
% degree_method (String)
%       'Relative' : Define polynomials in terms of degree with respect to
%                    x and y, so matrices of coefficients are rectangular.
%       'Total' :   Define polynomials in terms of total degree, so matrices
%                   of coefficients are square, where lower right triangle
%                   are all zero.
%
%       'Both' :    Combination of both of the above, and typically gives best
%                   results.
%
% % Example
%
% >> o_gcd_Bivariate_3Polys('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'None', 'None', 'Both')
%
% >> o_gcd_Bivariate_3Polys('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'Standard STLN','Standard Nonlinear APF','Relative')


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
SetGlobalVariables(problem_type, ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, degree_method)

% Print Parameters to screen
fprintf('INPUTS. \n')
fprintf('EXAMPLE NUMBER %s \n',ex_num)
fprintf('EMIN : %s \n',emin)
fprintf('EMAX : %s \n',emax)
fprintf('MEAN METHOD : %s \n', mean_method)
fprintf('PREPROCESSING : %s \n',num2str(bool_alpha_theta))
fprintf('LOW RANK METHOD : %s \n',low_rank_approx_method)
fprintf('APF METHOD : %s \n', apf_method)
fprintf('DEGREE METHOD : %s \n', degree_method)

% Get example polynomials
[fxy_exact, gxy_exact, hxy_exact,...
    uxy_exact,vxy_exact, wxy_exact,...
    dxy_exact,...
    m,m1,m2,...
    n,n1,n2,...
    o,o1,o2,...
    t_exact, t1_exact, t2_exact] = Examples_GCD_Bivariate_3Polys(ex_num);




DisplayDegreeStructure_3Polys();

% %
% %
% Add Noise

% Add noise to the coefficients of f and g
[fxy, ~] = AddVariableNoiseToPoly(fxy_exact, emin, emax);
[gxy, ~] = AddVariableNoiseToPoly(gxy_exact, emin, emax);
[hxy, ~] = AddVariableNoiseToPoly(hxy_exact, emin, emax);

% %
% % Get the GCD by zengs method
%[u,v,w] = o_gcd_zeng(fxy,gxy);


% Get GCD d(x,y) and quotient polynomials u(x,y) and v(x,y)
lowerLimit_t = 1;
upperLimit_t = min([m,n,o]);
limits_t = [lowerLimit_t,upperLimit_t];

lowerLimit_t1 = 0;
upperLimit_t1 = min([m1, n1, o1]);
limits_t1 = [lowerLimit_t1, upperLimit_t1];

lowerLimit_t2 = 0;
upperLimit_t2 = min([m2, n2, o2]);
limits_t2 = [lowerLimit_t2, upperLimit_t2];



switch SETTINGS.DEGREE_METHOD
    case 'Total'
        fxy_matrix_padd = zeros(m+1,m+1);
        gxy_matrix_padd = zeros(n+1,n+1);
        
        [r,c] = size(fxy);
        fxy_matrix_padd(1:r,1:c) = fxy;
        
        [r,c] = size(gxy);
        gxy_matrix_padd(1:r,1:c) = gxy;
        
        fxy = fxy_matrix_padd;
        gxy = gxy_matrix_padd;
        
    case 'Relative'
        
    case 'Both'
end

% Get the GCD by my method
[fxy_calc, gxy_calc, hxy_calc, dxy_calc, uxy_calc, vxy_calc, wxy_calc, t,t1,t2] = ...
    o_gcd_mymethod_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, limits_t, limits_t1, limits_t2);




% Compare exact d(x,y) and calculated d(x,y)
% PrintCoefficients(dxy_calc,dxy_exact,'d(x,y)');

switch SETTINGS.DEGREE_METHOD
    case 'Relative'
        
        error.dxy = GetDistanceBetweenPolynomials(dxy_exact, dxy_calc, 'd(x,y)');
        error.uxy = GetDistanceBetweenPolynomials(uxy_exact, uxy_calc, 'u(x,y)');
        error.vxy = GetDistanceBetweenPolynomials(vxy_exact, vxy_calc, 'v(x,y)');
        error.wxy = GetDistanceBetweenPolynomials(wxy_exact, wxy_calc, 'w(x,y)');
        
    case 'Total'
        
        % Get d(x,y) in a matrix in terms of total degree.
        dxy_exact_total = zeros(t_exact+1,t_exact+1);
        dxy_exact_total(1:t1_exact+1,1:t2_exact+1) = dxy_exact;
        
        uxy_exact_total = zeros(m-t_exact+1, m-t_exact+1);
        uxy_exact_total(1:m1-t1_exact+1,1:m2-t2_exact+1) = uxy_exact;
        
        vxy_exact_total = zeros(n-t_exact+1, n-t_exact+1);
        vxy_exact_total(1:n1-t1_exact+1,1:n2-t2_exact+1) = vxy_exact;
        
        wxy_exact_total = zeros(o-t_exact+1, o-t_exact+1);
        wxy_exact_total(1:o1-t1_exact+1,1:o2-t2_exact+1) = wxy_exact;
        
        error.dxy = GetDistanceBetweenPolynomials(dxy_exact_total, dxy_calc, 'd(x,y)');
        error.uxy = GetDistanceBetweenPolynomials(uxy_exact_total, uxy_calc, 'u(x,y)');
        error.vxy = GetDistanceBetweenPolynomials(vxy_exact_total, vxy_calc, 'v(x,y)');
        error.wxy = GetDistanceBetweenPolynomials(wxy_exact_total, wxy_calc, 'w(x,y)');
        
    case 'Both'

        error.dxy = GetDistanceBetweenPolynomials(dxy_exact,dxy_calc,'d(x,y)');
        error.uxy = GetDistanceBetweenPolynomials(uxy_exact,uxy_calc,'u(x,y)');
        error.vxy = GetDistanceBetweenPolynomials(vxy_exact,vxy_calc,'v(x,y)');
        error.wxy = GetDistanceBetweenPolynomials(wxy_exact,wxy_calc,'w(x,y)');
end



PrintToFile_GCD_3Polys(m, n, o, t, t1, t2, error)

% Given the two polynomials f(x,y) and g(x,y), Plot the explicit surfaces
% z = f(x,y) and z = g(x,y).
% switch SETTINGS.PLOT_GRAPHS
%     case 'y'
%         
%         surfaces{1} = fxy;
%         surfaces{2} = gxy;
%         surfaces{3} = hxy;
%         surfaces{4} = dxy_calc;
%         
%         PlotExplicitSurfaces(surfaces);
%         
%     case 'n'
%         
%     otherwise
%         error('error - plotgraphs is either y or n');
% end




end





function []= PrintToFile_GCD_3Polys(m, n, o, t, t1, t2, error)
%
% % Inputs
%
% m : (Int) Degree of polynomial f(x,y)
%
% n : (Int) Degree of polynomial g(x,y)
%
% o : (Int) Degree of polynomial h(x,y)
%
% t : (Int) Degree of polynomial d(x,y)
%
% t1 : (Int) Degree of polynomial d(x,y) with respect to x
%
% t2 : (Int) Degree of polynomial d(x,y) with respect to y
%
% error : Error.uxy Error.vxy Error.wxy Error.dxy


global SETTINGS

fullFileName = sprintf('Results/Results_o_gcd_3Polys%s.txt',datetime('today'));

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
       fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
        datetime('now'),... 
        SETTINGS.EX_NUM,...
        int2str(m),...
        int2str(n),...
        int2str(o),...
        int2str(t),...
        int2str(t1),...
        int2str(t2),...
        num2str(error.uxy),...
        num2str(error.vxy),...
        num2str(error.wxy),...
        num2str(error.dxy),...
        SETTINGS.MEAN_METHOD,...
        SETTINGS.BOOL_ALPHA_THETA,...
        SETTINGS.EMIN,...
        SETTINGS.EMAX,...
        SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
        int2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
        SETTINGS.APF_METHOD,...
        SETTINGS.DEGREE_METHOD); 
    end


    function WriteHeader()
        fprintf(fileID,'Date,EX_NUM,m,n,o,t,t1,t2,error_uxy,error_vxy,error_wxy,error_dxy,MEAN_METHOD,BOOL_ALPHA_THETA,EMIN,EMAX,LOW_RANK_APPROXIMATION_METHOD,ITERATIONS,APF_METHOD,DEGREE_METHOD\n');
    end








end
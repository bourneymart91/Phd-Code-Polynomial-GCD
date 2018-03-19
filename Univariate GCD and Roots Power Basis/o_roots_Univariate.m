function [] = o_roots_Univariate(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, deconvolution_method)
%
% O_ROOTS_UNIVARIATE(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
%
%
% Calculate the roots of a polynomial f(x) by a series of methods.
%
% % Inputs.
%
% ex_num : Example Number
%          Note a cusom example can be specified using the reg exp
%          'Custom:m=(\d+).low=(-?\d+).high=(-?\d+)'
%
% emin : Lower Noise Level
%
% emax : Upper Noise Level
%
% mean_method : Method of mean used to divide coefficients of Sylvester
%               Matrix as part of preprocessing.
%   'None' : No mean method is used
%   'Geometric Mean Matlab Method' : Divide by Geometric Mean
%
% bool_alpha_theta :
%   true : Include Preprocessing
%   false : Exclude Preprocessing
%
% low_rank_approx_method :
%   'None'
%   'Standard STLN'
%   'Standard SNTLN'
%   'Root Specific SNTLN'
%
% apf_method
%   'None'
%   'Standard APF Nonlinear' - Not Developed
%   'Standard APF Linear' - Not Developed
%
% % Example
% >> O_ROOTS_UNIVARIATE('1', 1e-10, 1e-12, 'None', false, 'None', 'None', 'Seperate')
% >> O_ROOTS_UNIVARIATE('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'None','None', 'Seperate')
% >> O_ROOTS_UNIVARIATE('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'Standard STLN','Standard APF Nonlinear', 'Seperate')
% >> O_ROOTS_UNIVARIATE('Custom:m=5 low=-1 high=2', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'Standard STLN', 'Standard APF Nonlinear', 'Seperate')

% Add Subfolders
restoredefaultpath
addpath(genpath(pwd));


% Initialise global variables
global SETTINGS

if ~isempty(SETTINGS)
else
    % If global variables are not yet set
    SETTINGS.DECONVOLUTION_METHOD = 'Batch';
    SETTINGS.ROOTS_UX = 'From ux';
end

% Set the problem type to be of type 'ROOTS'
problem_type = 'Roots';

% Set global input by user.
SetGlobalVariables_RootFinding(problem_type, ex_num, emin, emax, ...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, ...
    deconvolution_method)

% %
% %
% Get polynomial f(x)
fx_exact = Examples_Roots(ex_num);

% Add noise to coefficients of f(x)
fx = AddNoiseToPoly(fx_exact,emin);

% Create an array of methods used to compute polynomial factorisation
arr_factorisation_methods = ...
    {...
    'My Method',...
    'Musser Method'...
    'Matlab Method',...
    'Zeng Method'...
    };

% Get number of methods used
nMethods = length(arr_factorisation_methods);

% Initialise arrays
arr_root_multiplicity_matrix = cell(nMethods, 1);
arr_error = cell(nMethods, 1);
arr_time = cell(nMethods, 1);

for i = 1 : 1 : nMethods
    
    myTimer = tic();
    
    % Get method name
    method_name = arr_factorisation_methods{i};
    
    % Get roots by chosen method
    arr_root_multiplicity_matrix{i} = GetRoots(fx, method_name);
    PrintRoots(arr_root_multiplicity_matrix{i}, method_name);
    
    % Ger error
    arr_error{i} = GetRelativeError(arr_root_multiplicity_matrix{i}, fx_exact, method_name);
    
    % Get time taken
    arr_time{i} = toc(myTimer);
    
end

% Print to results file
PrintToFile(arr_factorisation_methods, arr_error, arr_time)

%
% Plot the graph real (x) imaginary (y) components of the nondistinct roots
% obtained by the root calculating methods.
if( SETTINGS.PLOT_GRAPHS)
    
    figure_name = sprintf('%s : Plot Calculated Roots', mfilename);
    figure('name', figure_name)
    hold on;

    for i = 1 : 1 : nMethods
        
        method_name = arr_factorisation_methods{i};
        
        root_multiplicity_matrix = arr_root_multiplicity_matrix{i};
        
        scatter( real(root_multiplicity_matrix(:,1)), imag(root_multiplicity_matrix(:,1)),'*','DisplayName', method_name);
        
    end
    
    
    
    xlabel('Real');
    ylabel('Imaginary');
    legend(gca,'show')
    ylim()
    str = sprintf('Plot of Calculated Roots of Polynomial f(y). \n componentwise noise = %g',emin);
    title(str);
    hold off
    
end


end

function [rel_err] = GetRelativeError(root_multiplicity_matrix, fx_exact, method)
% Get distance between the computed polynomial, and the exact polynomial,
% given the computed roots.
%
% Inputs.
%
% root_multiplicity_matrix : (Matrix)
%
% fx_exact : (Vector)
%
% method : (String)

try
    
    fx_computed = GetCoefficientsFromRoots(root_multiplicity_matrix);
    
    fprintf('Distance between f_exact and f_comp by %s : \n',method)
    
    % Get the relative error.
    rel_err = norm(fx_exact - fx_computed) ./ norm(fx_exact) ;
    
catch
    
    rel_err = 100000000000000;
    
end
end


function [] = PrintToFile(arr_Methods, arr_error, arr_time)
%
% % Inputs
%
% arr_Methods : (Array of String)
%
% arr_error : (Array of Float)


global SETTINGS

% Get number of methods used to compute the factorisation
nMethods = length(arr_Methods);

% Generate File name
fullFileName = sprintf('Results/Results_o_roots.txt');

% Check if file already exists
if exist(fullFileName, 'file')
    fileID = fopen(fullFileName,'a');
else
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
end

for i = 1 : 1 : nMethods
    
    method_name = arr_Methods{i};
    my_error = arr_error{i};
    my_time = arr_time{i};
    writeNewLine(method_name, my_error, my_time)
    
end


fclose(fileID);



    function writeNewLine(my_method_name, my_error, my_time)
        format_str = ...
            '%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s \n';
        
        fprintf(fileID,format_str,...
            datetime('now'),...
            SETTINGS.EX_NUM,...
            my_error,...
            my_time,...
            my_method_name,...
            SETTINGS.MEAN_METHOD,...
            num2str(SETTINGS.BOOL_ALPHA_THETA),...
            SETTINGS.EMIN,...
            SETTINGS.EMAX,...
            SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
            SETTINGS.DECONVOLUTION_METHOD_FX_HX,...
            SETTINGS.ROOTS_HX_COMPUTATION_METHOD...
            );
    end

    function WriteHeader()
        fprintf(fileID,'DATE, EX_NUM, ERROR, TIME, METHOD, MEAN_METHOD, PREPROC, EMIN, EMAX, LRA, DECONV METHOD, HX_METHOD \n');
    end

end


function [root_multiplicity_matrix] = GetRoots(fx,method_name)
%
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% method_name : (String)
%
% % Outputs
%
% root_multiplicity_matrix : (Matrix)


switch method_name
    
    case 'My Method'
        
        [root_multiplicity_matrix] = o_roots_mymethod(fx);
        
    case 'Musser Method'
        
        [root_multiplicity_matrix] = o_roots_Musser(fx);
        
    case 'Yun Method'
        
        [root_multiplicity_matrix] = o_roots_Yun(fx);
        
    case 'Matlab Method'
        
        [root_multiplicity_matrix] = o_roots_matlab(fx);
        
    case 'Zeng Method'
        
        [root_multiplicity_matrix] = o_roots_multroot(flipud(fx));
        
    otherwise
        
        error('Not a valid method')
        
end

end
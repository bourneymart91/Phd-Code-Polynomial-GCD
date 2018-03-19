function [] = o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method)
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
% >> O_ROOTS_UNIVARIATE('1', 1e-10, 1e-12, 'None', false, 'None','None')
% >> O_ROOTS_UNIVARIATE('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'None','None')
% >> O_ROOTS_UNIVARIATE('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'Standard STLN','Standard APF Nonlinear')
% >> O_ROOTS_UNIVARIATE('Custom:m=5 low=-1 high=2', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'Standard STLN','Standard APF Nonlinear')

% Add Subfolders
restoredefaultpath
addpath(...
    'Build Matrices',...
    'Deconvolution',...
    'Formatting',...
    'GCD Finding',...
    'Get Cofactor Coefficients',...
    'Get GCD Coefficients',...
    'Low Rank Approximation',...
    'Plotting',...
    'Preprocessing'...
    );

addpath(genpath('APF'));
addpath(genpath('Examples'));
addpath(genpath('Get GCD Degree'));
addpath(genpath('Root Finding'));
addpath(genpath('Low Rank Approximation'));

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
SetGlobalVariables(problem_type, ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method)

% %
% %
% Get polynomial f(x)
fx = Examples_Roots(ex_num);

% Add noise to coefficients of f(x)
fx = AddNoiseToPoly(fx,emin);

% %
% % MY METHOD
% %
% %

arr_factorisation_methods = {...
    'My Method',...
    'Musser Method'...
    'Yun Method',...
    'Matlab Method',...
    'Zeng Method'}


arr_root_mult_matrix{i}

for i = 1 : 1 : length(arr_factorisation_methods)
    
    % Get method name
    method_name = arr_factorisation_methods{i};
    
    % Get roots by chosen method
    [arr_root_mult_matrix{i}] = GetRoots(fx, method_name);
    
    GetRelativeError(arr_root_mult_matrix{i}
    
    
end





% %
% % MATLAB
% Get roots by matlab method

rel_err.MatlabMethod = GetRelativeError(root_mult_array_MatlabMethod,fx,'Matlab Method');
LineBreakLarge()


% %
% % MULTROOT
% Get roots by zheng method
computed_root_mult_array_multroot = o_roots_multroot(flipud(fx));
rel_err.MultrootMethod = GetRelativeError(computed_root_mult_array_multroot,fx,'MultRoot Method');
LineBreakLarge()


PrintToFile(rel_err,time)

%
% Plot the graph real (x) imaginary (y) components of the nondistinct roots
% obtained by the root calculating methods.
if( SETTINGS.PLOT_GRAPHS)
    
    figure_name = sprintf('%s : Plot Calculated Roots',mfilename);
    figure('name',figure_name)
    hold on;
    scatter( real(root_multiplicity_array_mymethod(:,1)), imag(root_multiplicity_array_mymethod(:,1)),'yellow','*','DisplayName','My Method');
    scatter( real(root_mult_array_MatlabMethod(:,1)), imag(root_mult_array_MatlabMethod(:,1)),'red','DisplayName','Matlab Roots');
    scatter( real(computed_root_mult_array_multroot(:,1)), imag(computed_root_mult_array_multroot(:,1)),'green','s','filled','DisplayName','MultRoots');
    
    
    xlabel('Real');
    ylabel('Imaginary');
    legend(gca,'show')
    ylim()
    str = sprintf('Plot of Calculated Roots of Polynomial f(y). \n componentwise noise = %g',emin);
    title(str);
    hold off
    
end


end

function [rel_err] = GetRelativeError(computed_root_mult_array,fx_exact,method)
% Get distance between the computed polynomial, and the exact polynomial,
% given the computed roots.
%
% Inputs.
%
% computed_root_mult_array
%
% method : (String)

try
    fx_computed = GetCoefficients(computed_root_mult_array);
    
    fprintf('Distance between f_exact and f_comp by %s : \n',method)
    
    % Get the relative error.
    rel_err = norm(fx_exact - fx_computed) ./ norm(fx_exact) ;
    
    display(rel_err);
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


nMethods = length(arr_Methods);

[y, m, d] = ymd(datetime('now'));

fullFileName = sprintf('Results/Results_o_roots.txt_%s-%s-%s',y,m,d);


if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    
    format_str = ...
        '%s, %s, %s, %s, %s, %s, %s, %s, %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s \n';
    
    
    for i = 1:1:nMethods

        my_method_name = arr_Methods{i};
        my_error = arr_error{i};
        my_time = arr_time{i};
        
        
        fprintf(fileID,format_str,...
            datetime('now'),...
            SETTINGS.EX_NUM,...
            my_error,...
            my_time,...
            my_method_name,...
            SETTINGS.MEAN_METHOD,...
            SETTINGS.BOOL_ALPHA_THETA,...
            SETTINGS.EMIN,...
            SETTINGS.EMAX,...
            SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
            SETTINGS.DECONVOLUTION_METHOD_FX_HX,...
            SETTINGS.ROOTS_HX_COMPUTATION_METHOD...
            );

    end
    
    fclose(fileID);
    
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
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
        
    otherwise
        
        error('Not a valid method')
        
end

end
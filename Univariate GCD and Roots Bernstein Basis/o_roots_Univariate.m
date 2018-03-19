function [epsilon_fi, epsilon_hi, epsilon_wi] = o_roots_Univariate(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_variant, rank_revealing_metric, deconvolution_method_hx,...
    deconvolution_method_wx, deconvolution_preproc)
% O_ROOTS_UNIVARIATE(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
%   low_rank_approx_method, apf_method, sylvester_matrix_variant, ...
%   rank_revealing_metric, deconvolution_method_hx, deconvolution_method_wx)
%
% Given an example number, and a set of input parameters, calculate the
% roots r_{i} of the polynomial f(x) and the corresponding multiplicities
% m_{i}.
%
% % Inputs
%
% ex_num : (String) Example Number
%
% emin : (Float) Noise/Signal maximum threshold (minimum)
%
% emax : (Float) Noise/Signal maximum threshold (maximum)
%
% mean_method : (string) method used to compute the mean of entries in C_{n-k}(f)
%               and C_{m-k}(g)
%           'None' - No mean
%           'Geometric Mean Matlab Method'
%           'Geometric Mean My Method'
%
%
% bool_alpha_theta : (Boolean)
%           true
%           false
%
% low_rank_approx_method : (string)
%           'None'
%           'Standard STLN'
%           'Standard SNTLN'
%           'Root Specific SNTLN'
%
% apf_method ('string')
%           'None'
%           'Standard APF NonLinear'
%           'Standard APF Linear'
%
% sylvester_matrix_variant
%   * 'T'
%   * 'DT'
%   * 'DTQ'
%   * 'TQ'
%   * 'DTQ Denominator Removed'
%   * 'DTQ Rearranged'
%
% rank_revealing_metric
%   * Singular Values
%   * Max/Min Singular Values
%   * R1 Row Norms
%   * R1 Row Diagonals
%   * Residuals
%
% deconvolution_method_hx
%   * Separate
%   * Batch
%   * Batch With STLN
%   * Batch Constrained
%   * Batch Constrained With STLN
%
%
% deconvolution_method_wx
%   * Separate
%   * Batch
%
%
% % Example
% >> O_ROOTS_UNIVARIATE('1', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'None', 'None', 'DTQ', 'Minimum Singular Values', 'Batch Constrained', 'Batch', true)



% Add subfolders
restoredefaultpath
addpath(genpath(pwd));




% Set the global variables
global SETTINGS

SetGlobalVariables_Roots(ex_num, emin, emax,...
    mean_method, bool_alpha_theta, low_rank_approx_method, ...
    apf_method, sylvester_matrix_variant, rank_revealing_metric, ...
    deconvolution_method_hx, deconvolution_method_wx, deconvolution_preproc)



% Check that max and min signal to noise ratio are the correct way around.
% If not, rearrange min and max.
if emin > emax
    
    fprintf('Minimum noise greater than maximum noise \n Swapping values...\n')
    emin_wrong = emin;
    emax_wrong = emax;
    emin = emax_wrong;
    emax = emin_wrong;
    fprintf('')
    
end


% Print Settings to console
PrintGlobalVariables();

% Get the polynomial f(x) as a column vector of coefficients.
[fx_exact, fx_factor_multiplicity_matrix] = Examples_Roots(ex_num);

% Get Array of polynomials f_{i}(x)
arr_fx_exact = GetArray_fx(fx_factor_multiplicity_matrix);

% Get Array of polynomials h_{i}(x)
arr_hx_exact = GetArray_hx(fx_factor_multiplicity_matrix);

% Get Array of polynomials w_{i}(x)
arr_wx_exact = GetArray_wx(fx_factor_multiplicity_matrix);

% Add Noise to coefficients of exact polynomial f_exact, to obtain noisy
% polynomial fx.
fx = AddVariableNoiseToPoly(fx_exact, emin, emax);




% Root Finding Methods
%   My Method
%   Musser Method
%   Yun Method
%   Matlab Method

arr_RootFindingMethod = {...
    'My Method', ...
    'Matlab Method', ...
    'Zeng Method',...
    %'Musser Method', ...
    %'Bisection Method', ...
    %'Subdivision Method', ...
    %'Bezier Clipping Method'...
    };
%arr_RootFindingMethod = {'Yun Method'};


% Initialise some arrays
nMethods = length(arr_RootFindingMethod);
arrRootMultiplicity = cell(nMethods,1);

arrForwardErrors = cell(nMethods,1);
arrBackwardErrors = cell(nMethods,1);



% For each method, compute the roots and multiplicities



for i = 1 : 1 : nMethods
    
    % Get name of root finding Method
    method_name = arr_RootFindingMethod{i};
    
    LineBreakLarge();
    fprintf(['Method : ' method_name '\n']);
    LineBreakLarge();
    
    
    
    
    if strcmp(method_name, 'My Method')
        
        % Get roots by my method
        [arrRootMultiplicity{i}, arr_fx_comp_myMethod, ...
            arr_hx_comp_myMethod, arr_wx_comp_myMethod] = o_roots_mymethod(fx);
        
        
      
        
        % Compare exact polynomials f_{i}(x) and computed polynomials
        % f_{i}(x)
        try
            vErrors_arr_fx = GetPolynomialArrayErrors(arr_fx_exact, arr_fx_comp_myMethod);
            %PlotErrors_fi(vErrors_arr_fx, low_rank_approx_method);
            
        catch
            vErrors_arr_fx = 1000;
        end
        
        
        epsilon_fi = (mean(vErrors_arr_fx));
        fprintf('Average Error in f_{i}(x) : %e \n', epsilon_fi);
        
        
        
        % Compare exact polynomials h_{i}(x) and computed polynomials
        % h_{i}(x)
        
        try
            vErrors_arr_hx = GetPolynomialArrayErrors(arr_hx_exact, arr_hx_comp_myMethod);
            %PlotErrors_hi(vErrors_arr_hx, deconvolution_method_hx)
            
        catch
            
            vErrors_arr_hx = 1000;
        end
        
        epsilon_hi = mean((vErrors_arr_hx));
        fprintf('Average Error in h_{i}(x) : %e \n', epsilon_hi);
        
        
        
        
        
        % Compare exact polynomials w_{i}(x) with computed polynomials
        % w_{i}(x)
        
        try
            
            if length(arr_wx_exact) == length(arr_wx_comp_myMethod)
            
                vErrors_arr_wx = GetPolynomialArrayErrors(arr_wx_exact, arr_wx_comp_myMethod);
                %PlotErrors_wi(vErrors_arr_wx, deconvolution_method_wx)
            
            else
                
                vErrors_arr_wx = 1000 * ones(length(arr_wx_exact),1);
                
            end
            
        catch
            
            vErrors_arr_wx = 1000;
        end
        
        epsilon_wi = mean((vErrors_arr_wx(vErrors_arr_wx~=0)));
        fprintf('Average Error in w_{i}(x) : %e \n', epsilon_wi);
        
        
        
        %try
            
            % Plot errors in h_{i}(x)
            
            
            figure_name = sprintf('Errors in f_{i},h_{i}, w_{i}(x) ');
          
            line_name_h = sprintf('$Deconvolution h_{i}(x) Method : %s - %s$', deconvolution_method_hx, num2str(deconvolution_preproc));
            line_name_w = sprintf('$Deconvolution w_{i}(x) Method : %s$', deconvolution_method_wx);
            line_name_f = sprintf('$Low Rank Approx : %s$', low_rank_approx_method);
            
            figure('Name', figure_name)
            
            hold on
            
            plot(log10(vErrors_arr_fx), '-o', 'DisplayName', line_name_f, 'Color', 'green', 'MarkerFaceColor', 'green')
            plot(log10(vErrors_arr_hx), '-o', 'DisplayName', line_name_h, 'Color', 'blue', 'MarkerFaceColor', 'blue')
            plot(log10(vErrors_arr_wx), '-o', 'DisplayName', line_name_w, 'Color', 'red', 'MarkerFaceColor', 'red')
            
            % Labels and Legends
            xlabel('$i$','Interpreter','latex')
            ylabel('$\log_{10} \left( \epsilon \right)$','Interpreter', 'latex')

            l = legend(gca,'show');
            set(l,{'Location', 'Interpreter'},{'southeast', 'latex'});
            
            grid on
            box on
            hold off
            
            
        %catch
            
        %end
        
        
        
    else
        
        % Get roots and multiplicities
        arrRootMultiplicity{i} = GetRootsAndMultiplicities(fx, method_name);
        
    end
    % Get Forward and backward error
    try
        arrBackwardErrors{i} = GetBackwardErrorMeasure(arrRootMultiplicity{i}, fx_exact);
        arrForwardErrors{i} = GetForwardErrorMeasure(arrRootMultiplicity{i}, fx_factor_multiplicity_matrix);
        
    catch
        arrForwardErrors{i} = 1000;
    end
    
    
    
    LineBreakLarge()
    
end


% Print Roots to command line
LineBreakLarge()
LineBreakLarge()

for i = 1 : 1 : nMethods
    
    method_name = arr_RootFindingMethod{i};
    root_mult_matrix = arrRootMultiplicity{i};
    backward_error = arrBackwardErrors{i};
    forward_error = arrForwardErrors{i};
    
    
    PrintoutRoots(method_name, root_mult_matrix);
    PrintErrors(method_name, backward_error, forward_error);
    
    
    
end




%
% Plot the graph real (x) imaginary (y) components of the nondistinct roots
% obtained by the root calculating methods.
if (SETTINGS.PLOT_GRAPHS_ROOTS)
    
    figure_name = sprintf('%s : Plot Calculated Roots',mfilename);
    figure('name',figure_name)
    hold on;
    
    
    for i = 1:1:length(arr_RootFindingMethod)
        
        % Get method name
        methodName = arr_RootFindingMethod{i};
        
        % Get matrix of roots and corresponding multiplicities
        mat_Root_Mult = arrRootMultiplicity{i};
        
        try
            % Get Vector of real part of each root
            vRealPart = real(mat_Root_Mult(:,1));
            
            % Get Vector of imag part of each root
            vImagPart = imag(mat_Root_Mult(:,1));
            
            % Get scatter data
            scatter( vRealPart, vImagPart, 'DisplayName', methodName);
            
        catch
            fprintf('Roots found by %s method can not be plotted. \n', methodName);
        end
        
    end
    
    
    
    grid on
    xlabel('Real');
    ylabel('Imaginary');
    legend(gca,'show');
    ylim();
    str = sprintf('Plot of Calculated Roots of Polynomial f(y). \n componentwise noise = %g',emin);
    title(str);
    hold off
    
end


PrintRootsToFile(arr_RootFindingMethod, arrBackwardErrors, arrForwardErrors, vErrors_arr_fx, vErrors_arr_hx, vErrors_arr_wx);


end

function PrintErrors(method_name, backward_error, forward_error)

fprintf('Method Name : %s \n', method_name)
fprintf('Foward Error : %e \n', forward_error);
fprintf('Backward Error : %e \n', backward_error);

end








function [] = PlotErrors_wi(vErrors_arr_wx, deconvolution_method_wx)
% Plot errors in h_{i}(x)
figure_name = sprintf('Errors in w_{i}(x) %s', deconvolution_method_wx);
line_name_w = sprintf('Method : %s', deconvolution_method_wx);


figure('Name', figure_name)
hold on
plot(log10(vErrors_arr_wx),'-s','DisplayName', line_name_w)
xlabel('$w_{i}(x)$','Interpreter', 'latex')
ylabel('$\log_{10} \left( \epsilon w_{i}(x) \right) $', 'Interpreter', 'latex')
hold off
end


function [] = PlotErrors_hi(vErrors_arr_hx, deconvolution_method_hx)


% Plot errors in h_{i}(x)
figure_name = sprintf('Errors in h_{i}(x) %s', deconvolution_method_hx);
line_name_h = sprintf('Method : %s', deconvolution_method_hx);


figure('Name', figure_name)
hold on
plot(log10(vErrors_arr_hx),'-s','DisplayName', line_name_h)
xlabel('$h_{i}(x)$','Interpreter', 'latex')
ylabel('$\log_{10} \left( \epsilon h_{i}(x) \right)$', 'Interpreter', 'latex')
hold off


end


function [] = PlotErrors_fi(vErrors_arr_fx, low_rank_approx_method)

% Plot errors
figure_name = sprintf('Errors in f_{i}(x) - LRA : %s', low_rank_approx_method);
line_name_f = sprintf('%s', low_rank_approx_method);

figure('Name', figure_name)
    hold on
    plot(log10(vErrors_arr_fx),'-s','DisplayName',line_name_f)
    xlabel('$f_{i}(x)$', 'Interpreter', 'latex')
    ylabel('$\log_{10} \left( \epsilon f_{i}(x) \right)$', 'Interpreter', 'latex')
    hold off
end
function [] = o_Deconvolution(ex_num, emin, bool_preproc)
% Test the different methods of deconvolving the polynomials f_{i}(x), to
% form the set of polynomials h_{i} where h_{i} = f{i}/f_{i+1}
%
% % Inputs
%
% ex_num : (String) Example number (String)
%
% noise : (Float) noise level
%
% bool_preproc : (Boolean) Bool determining whether to include preprocessing
%
%
%
% % Outputs
%
%
% Results are printed to a dat file
%
%
% % Example
%
% >> o_Deconvolution('1', 1e-12, true)


% Set global settings
global SETTINGS
SETTINGS.PLOT_GRAPHS = true;
SETTINGS.PLOT_GRAPHS_DECONVOLUTION_LRA = false;

SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-12;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 50;

SETTINGS.BOOL_LOG = false;
SETTINGS.PREPROC_DECONVOLUTIONS = bool_preproc;

SETTINGS.SEED = 1024;


restoredefaultpath();
addpath(genpath('../Examples'));

% Add that folder plus all subfolders to the path.
addpath(genpath(pwd));


syms x y;


% Get the array containing symbolic factors and their corresponding
% multiplicity in f_{0}(x)
[factor_multiplicity_matrix] = Deconvolution_Examples_Univariate(ex_num);

% Get array of polynomials f_{i}(x)
arr_fx_exact = GetArray_fx(factor_multiplicity_matrix);

% Get array of polynomials h_{i}(x)
arr_hx_exact = GetArray_hx(factor_multiplicity_matrix);

% Get number of polynomials in array f(x)
nPolynomials_arr_fx = length(arr_fx_exact);

% Get degree of polynomials h_{i}(x)
vDegree_arr_hx = GetDegree_Array(arr_hx_exact);

% Get the degree structure of the polynomials w_{i} where w_{i} =
% h_{i-1}/h_{i}
vDegree_arr_wx = diff([vDegree_arr_hx; 0]);

% Get the multiplicities of the factors of f(x)
vMultiplicities = find(vDegree_arr_wx~=0);






% Add noise to the coefficients of f_{i}(x)

% Initialise a cell array to store noisy polynomials f_{i}(x)
arr_fx_noisy = cell(nPolynomials_arr_fx, 1);

% Get noisy polynomials f_{i}(x)
for i = 1 : 1 : nPolynomials_arr_fx
    
    arr_fx_noisy{i,1} = AddNoiseToPoly(arr_fx_exact{i},emin);
    
end





% Define an array of deconvolution methods to be used
arr_DeconvolutionMethod = {...
    'Separate' ...
    'Batch', ...
    'Batch With STLN',...
    'Batch Constrained',...
    'Batch Constrained With STLN'...
    };


nMethods = length(arr_DeconvolutionMethod);

% Testing deconvolution
LineBreakLarge();
arr_hx = cell(nMethods,1);
arr_Error = cell(nMethods,1);

for i = 1 : 1 : nMethods
    
    % Get deconvolution method
    method_name = arr_DeconvolutionMethod{i};
    
    switch method_name
        
        case 'Separate'
            arr_hx{i,1} = Deconvolve_Separate(arr_fx_noisy);
            
        case 'Batch'
            arr_hx{i,1} = Deconvolve_Batch(arr_fx_noisy);
            
        case 'Batch With STLN'
            arr_hx{i,1} = Deconvolve_Batch_With_STLN(arr_fx_noisy);
            
        case 'Batch Constrained'
            arr_hx{i,1} = Deconvolve_Batch_Constrained(arr_fx_noisy, vMultiplicities);
            
        case 'Batch Constrained With STLN'
            arr_hx{i,1} = Deconvolve_Batch_Constrained_With_STLN(arr_fx_noisy, vMultiplicities);
            
        otherwise
            error('err')
            
    end
    
    fprintf([mfilename ' : ' sprintf('%s \n', method_name )]);
    
    arr_Error{i,1} = GetPolynomialArrayErrors(arr_hx{i,1}, arr_hx_exact);
    
end


PlotErrors(arr_Error, arr_DeconvolutionMethod)



%--------------------------------------------------------------------------
% Console writing

for i = 1:1:nMethods
    
    methodName = arr_DeconvolutionMethod{i};
    vError = arr_Error{i};
    display([mfilename ' : ' sprintf('Error %s : %e', methodName, mean(vError))]);
    
end



% Initialise array to store error for each method
arr_ErrorNorm = cell(nMethods,1);

for i = 1 : 1 : nMethods
    arr_ErrorNorm{i} = norm(arr_Error{i});
end


PrintToResultsFile(ex_num, bool_preproc, emin, arr_DeconvolutionMethod, arr_ErrorNorm);

end

function [] = PrintToResultsFile(ex_num, bool_preproc, noise, arr_DeconvolutionMethod, arr_ErrorNorm)
%
% % Inputs
%
% ex_num : (String)
%
% bool_preproc : (Boolean)
%
% noise : (float)
%
% arr_DeconvolutionMethod : (Array of Strings)
%
% arr_ErrorNorm : (Array of floats)


fullFileName = sprintf('Results/Results_o_deconvolutions.dat');
nMethods = length(arr_DeconvolutionMethod);

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    
    
    
    for i = 1:1:nMethods
        
        method_name = arr_DeconvolutionMethod{i};
        error_norm = arr_ErrorNorm{i};
        
        
        WriteNewLine(method_name, error_norm)
        
        
    end
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    
    WriteHeader()
    for i = 1 : 1 : nMethods
        
        method_name = arr_DeconvolutionMethod{i};
        error_norm = arr_ErrorNorm{i};
        WriteNewLine(method_name, error_norm);
        
    end
    fclose(fileID);
    
end

    function WriteNewLine(method_name, error_norm)
        
        %
        fprintf(fileID,'%s,%s,%s,%s,%s,%s \n',...
            datetime('now'),...
            ex_num,...
            num2str(bool_preproc),...
            num2str(noise),...
            method_name,...
            num2str(error_norm)...
            );
        
    end

    function WriteHeader()
        
        fprintf(fileID,'DATE, EX_NUM, BOOL_PREPROC, NOISE, method_name, error_norm \n');
        
    end


end





function PlotErrors(arr_Error, arr_DeconvolutionMethod)

% Plotting


global SETTINGS


nMethods = length(arr_DeconvolutionMethod);

if (SETTINGS.PLOT_GRAPHS)
    
    figure_name = sprintf([mfilename ':' 'Deconvolution Methods Error']);
    figure('name',figure_name)
    hold on
    
    
    for i = 1 : 1 : nMethods
        
        methodName = arr_DeconvolutionMethod{i};
        vError = arr_Error{i,1};
        plot(log10(vError), '-o', 'DisplayName', methodName)
        
    end
    
    nPolynomials_hx = length(vError);
    
    l = legend(gca,'show');
    set(l,'Interpreter','latex')
    
    
    
    % Resizing and Repositioning
    myplot = gca;
    myval_side = 0.10;
    myval_base = 0.08;
    set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])

    grid on
    box on
    
    xlim([1 nPolynomials_hx]);
    xlabel('i : Index of $f_{i}(x)$','Interpreter','latex','FontSize', 20)
    ylabel('$\log_{10}\left( \epsilon \right)$','Interpreter','latex','FontSize', 20)
    hold off
    
end
end



function [arr_fx] = GetArray_fx(fx_factor_multiplicity_matrix)
%
%
% fx_factor_multiplicity_matrix : Matrix containing the symbolic factors of
% f(x) and the multiplicity of the factors.
%
% % Outputs
%
% arr_fx : (Array of Vectors) Each vector contains the coefficients of a
% polynomial f_{i}(x) in the sequence generated in Gauss method of
% factorisation.

% Get symbolic factors
vFactors = (fx_factor_multiplicity_matrix(:,1));

% Get multiplicity of factors
vMultiplicity_f0 = double(fx_factor_multiplicity_matrix(:,2));

% Get max multiplicty
maxMultiplicity = max(vMultiplicity_f0);

% Set number of polynomials in array
nPolys_arr_fx = maxMultiplicity + 1;

% Initialise array
arrMultiplicities = cell(nPolys_arr_fx, 1);
arr_fx = cell(nPolys_arr_fx, 1);

for i = 1 : 1 : nPolys_arr_fx
    
    % Get multiplicity of factors in f_{i}(x)
    vMultiplicity_fi = vMultiplicity_f0 - (i-1);
    vMultiplicity_fi(vMultiplicity_fi < 0) = 0;
    
    arrMultiplicities{i} = vMultiplicity_fi;
    
    % Get coefficients of f_{i}(x)
    arr_fx{i} = BuildPolyFromRootsSymbolic([vFactors vMultiplicity_fi]);
    
end


end




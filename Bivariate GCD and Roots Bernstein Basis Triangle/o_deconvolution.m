function [] = o_deconvolution(ex_num, emin, bool_preproc)
% Perform tests on deconvolution methods
%
%
% % Inputs
%
% ex_num : (String) Example number (String)
%
% bool_preproc : (Boolean) whether preprocessing is included.
%
% emin : (Float) Noise level
%
%
% % Outputs.
%
%
%
% % Example
%
% >> o_deconvolution('1', true, 1e-10);

% add location for Bivariate_Deconvolution_Examples() function
addpath(genpath('../Examples'));

% Add that folder plus all subfolders to the path.
addpath(genpath(pwd));

% Add symbolic variables
syms x y;

global SETTINGS
SETTINGS.SEED = 1024;
SETTINGS.PREPROC_DECONVOLUTIONS = bool_preproc;
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-14;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 50;
SETTINGS.PLOT_GRAPHS = true;

% Get a matrix containing the symbolic factors of a polynomial and a vector
% containing the multiplicity of each factor.
[factor_mult_arr] = Deconvolution_Examples_Bivariate(ex_num);

vFactor = factor_mult_arr(:,1);

vMult = factor_mult_arr(:,2);

% vMult is symbolic so make it numeric
vMult = double(vMult);

% Get the highest multiplicity of any factor in f_{0}(x,y)
highest_multiplicity = max(vMult);


% %
% %
% Generate the polynomials f_{0},...,f_{m}(x)

% Initialise array to store polynomials f_{i}(x,y)
arr_symbolic_fxy = cell(highest_multiplicity + 1, 1);

% Initialise a vector to store the degree of each polynomial f_{i}(x,y)
vDegree_arr_fxy = zeros(highest_multiplicity + 1, 1);

% For each polynomial f_{i}(x,y)
for i = 0 : 1 : highest_multiplicity
    
    % Get the multiplicities of the roots of f_{i}
    mults = ((vMult - i) + abs(vMult - i)) ./ 2;
    
    % Get the symbolic polynomial f_{i}
    arr_symbolic_fxy{i + 1} = prod(vFactor.^(mults));
    
    % Get the degree of polynomial f_{i}(x)
    vDegree_arr_fxy(i + 1) = double(feval(symengine, 'degree', (arr_symbolic_fxy{i+1})));
    
end

% Display first polynomial
display(arr_symbolic_fxy{1})

%
%
%
%

% Get the degree structure of the polynomials h_{i} where h_{i} =
% f_{i-1}/f_{i}.
vDegree_arr_hxy = -1 .* diff(vDegree_arr_fxy);

% Get the sequence of polynomials h_{i}(x) in symbolic form
arr_symbolic_hxy = cell(length(arr_symbolic_fxy) - 1, 1);

for i = 1 : 1 : length(arr_symbolic_fxy) - 1
    
    arr_symbolic_hxy{i} = arr_symbolic_fxy{i} / arr_symbolic_fxy{i + 1};
    
end


% %
% %
% Get coefficients vectors of f_{i}(x) and h_{i}(x)
nPolynomials_arr_fxy = length(arr_symbolic_fxy);
nPolynomials_arr_hxy = length(arr_symbolic_fxy) - 1;

arr_fxy_power = cell(nPolynomials_arr_fxy, 1);
arr_fxy_bernstein = cell(nPolynomials_arr_fxy, 1);

arr_hxy_power_exact = cell(nPolynomials_arr_hxy, 1);
arr_hxy_bernstein_exact = cell(nPolynomials_arr_hxy, 1);

for i = 1:1:nPolynomials_arr_fxy
    
    
    arr_fxy_power{i, 1} = double(rot90(coeffs(arr_symbolic_fxy{i}, [x,y], 'All'), 2));
    
    arr_fxy_bernstein{i, 1} = PowerToBernstein(arr_fxy_power{i, 1}, vDegree_arr_fxy(i));
    
    if i <= nPolynomials_arr_hxy
        
        arr_hxy_power_exact{i, 1} = double(rot90(coeffs(arr_symbolic_hxy{i}, [x,y], 'All'), 2));
        
        arr_hxy_bernstein_exact{i, 1} = PowerToBernstein(arr_hxy_power_exact{i,1}, vDegree_arr_hxy(i));
        
    else
        arr_fxy_power{i, 1} = 1;
    end
    
end

% %
% %
% %
% Add noise to the coefficients of f_{i}(x)
arr_fxy_bernstein_noisy = cell(nPolynomials_arr_fxy,1);

for i = 1 : 1 : nPolynomials_arr_fxy
    
    arr_fxy_bernstein_noisy{i,1} = AddNoiseToPoly(arr_fxy_bernstein{i}, emin);
    
end

% Define an array of deconvolution methods to be used
arr_DeconvolutionMethod = {...
    'Separate',...
    'Batch',...
    'Batch With STLN',...
    'Batch Constrained',...
    'Batch Constrained With STLN'};


nMethods = length(arr_DeconvolutionMethod);

% Testing deconvolution
LineBreakLarge();
arr_hxy = cell(nMethods,1);
arr_Error = cell(nMethods,1);

% For each deconvolution method
for i = 1 : 1 : nMethods
    
    % Get name of deconvolution method
    method_name = arr_DeconvolutionMethod{i};
    
    switch method_name
        
        case 'Separate'
            arr_hxy{i,1} = Deconvolve_Bivariate_Separate(arr_fxy_bernstein_noisy);
            
        case 'Batch'
            arr_hxy{i,1} = Deconvolve_Bivariate_Batch(arr_fxy_bernstein_noisy, vDegree_arr_fxy);
            
        case 'Batch With STLN'
            arr_hxy{i,1} = Deconvolve_Bivariate_Batch_With_STLN(arr_fxy_bernstein_noisy, vDegree_arr_fxy);
            
        case 'Batch Constrained'
            arr_hxy{i,1} = Deconvolve_Bivariate_Batch_Constrained(arr_fxy_bernstein_noisy, vDegree_arr_fxy);
            
        case 'Batch Constrained With STLN'
            arr_hxy{i,1} = Deconvolve_Bivariate_Batch_Constrained_With_STLN(arr_fxy_bernstein_noisy, vDegree_arr_fxy);
            
        otherwise
            error('err')
            
    end
    
    fprintf([mfilename ' : ' sprintf('%s \n', method_name )]);
    
    arr_Error{i, 1} = GetErrors(arr_hxy{i}, arr_hxy_bernstein_exact);
    
end





% Plotting results

if (SETTINGS.PLOT_GRAPHS_DECONVOLUTION)
    
    figure_name = sprintf([mfilename ':' 'Deconvolution Methods Error']);
    figure('name',figure_name)
    hold on
    
    
    arr_styles= {'-s', '-o', '-+', '-*', '-square'};
    
    for i = 1 : 1 : nMethods
        
        methodName = arr_DeconvolutionMethod{i};
        vError = arr_Error{i};
        style = arr_styles{i};
        plot(log10(vError), style, 'DisplayName', methodName, 'LineWidth',2)
        
        
    end
    
    % Figure Layout
    
    grid on
    box on
    
    l = legend(gca,'show');
    set(l,'Interpreter','latex');
    set(l,'FontSize',16);
    
    
    myval_side = 0.10;
    myval_base = 0.08;
    myplot = gca;
    set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])
    set(gcf, 'Position', [100, 100, 710, 650])
    
    % Figure Labels
    xlabel('i : Index of factor $f_{i}(x)$', 'Interpreter', 'latex','FontSize',20);
    ylabel('$\log_{10} (\epsilon)$','Interpreter', 'latex', 'FontSize', 20);
    hold off
    
end


%--------------------------------------------------------------------------
% Console writing

for i = 1 : 1 : nMethods
    
    methodName = arr_DeconvolutionMethod{i};
    vError = arr_Error{i};
    display([mfilename ' : ' sprintf('Error %s : %e', methodName, mean((vError)))]);
    
end



% Initialise array to store error for each method
arr_ErrorNorm = cell(nMethods,1);

for i = 1 : 1 : nMethods
    arr_ErrorNorm{i} = norm(arr_Error{i});
end


PrintToResultsFile(ex_num, bool_preproc, emin, arr_DeconvolutionMethod, arr_ErrorNorm);


end

function [] = PrintToResultsFile(ex_num, bool_preproc, noise, arr_DeconvolutionMethod, arr_ErrorNorm)
% Print results of each deconvolution method to 
% 
% % Inputs
%
% ex_num : (String)
%
% bool_preproc : (Boolean)
%
% noise : (Float) 
%
% arr_DeconvolutionMethod : (Array of Strings)
%
% arr_ErrorNorm : (Array of Integers)
%



% Get number of deconvolution methods considered
nMethods = length(arr_DeconvolutionMethod);

fullFileName = sprintf('Deconvolution/Results/Results_o_deconvolutions%s.txt',datetime('today'));

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    
    for i = 1 : 1 : nMethods
        
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
        WriteNewLine(method_name, error_norm)
        
    end
    
    fclose(fileID);
    
end

    function WriteNewLine(method_name, error_norm)
        %
        % % Inputs
        %
        % method_name : (String)
        %
        % error_norm : (Float)
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
        
        fprintf(fileID,'DATE, EX_NUM, BOOL_PREPROC, NOISE, Method, Error \n');
        
    end


end


function [vError] = GetErrors(arr_fxy_exact, arr_hxy_approx)
%
% % Inputs
%
% arr_fxy_exact : (Array of Matrices)
%
% arr_hxy : (Array of Matrices)


vError = zeros(length(arr_fxy_exact),1);

% Compare each computed h{i} with actual h_{i}
for i = 1:1:length(arr_fxy_exact)
    
    
    hxy_approx = arr_hxy_approx{i};
    fxy_exact = arr_fxy_exact{i};
    
    
    vError(i) = GetMatrixDistance(fxy_exact, hxy_approx);
    
    
end

end




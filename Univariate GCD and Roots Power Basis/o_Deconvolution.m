function [] = o_Deconvolution(ex_num, emin, bool_preproc)
% O_DECONVOLUTION Test the different methods of deconvolving the set of
% polynomials f_{i}(x), to form the set of polynomials h_{i}
% where h_{i} = f{i}/f_{i+1}
%
% % Inputs
%
% ex_num : (String) Example number (String)
%
% emin : (Float) Lower noise level
%
% emax : (Float) Upper noise level
%
% % Outputs.
%
% Outputs are printed to file
%
% Example
%
% >> o_Deconvolution('1',1e-10, true)


% Add path for examples
restoredefaultpath;

% Add that folder plus all subfolders to the path.
addpath(genpath(pwd));
addpath(genpath('../Examples'));


% Set settings pertaining to this test
global SETTINGS
SETTINGS.PLOT_GRAPHS = true;
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-20;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 20;
SETTINGS.SEED = 1024;
SETTINGS.PREPROC_DECONVOLUTIONS = bool_preproc;



% % Get the factor array and multiplicity vector
[factor_mult_arr_f] = Deconvolution_Examples_Univariate(ex_num);

% Get the symbolic factors of f(x)
vFactors = factor_mult_arr_f(:,1);

% Get the multiplicities of the factors of f(x)
vMultiplicity = double(factor_mult_arr_f(:,2));

% Get highest power of any of the factors of f(x)
highest_pwr = max(vMultiplicity);

% Generate polynomials f_{0}(x) ,..., f_{m}(x) = 1. Where each f_{i+1}(x) is
% the f_{i+1} = GCD(f_{i},f'_{i}).
arr_sym_fx = cell(highest_pwr+1, 1);
vDegree_fx = zeros(highest_pwr+1, 1);

% Get the number of polynomials in the array f_{i}(x)
nPolys_arr_fx = highest_pwr + 1;

for i = 1 : 1 : nPolys_arr_fx
    
    % Get the multiplicities of the factors of f_{i+1}(x). The
    % multiplicities of each factor are one less than the previous
    % f_{i}(x).
    vMultiplicities = ((vMultiplicity - (i-1)) + abs(vMultiplicity - (i-1))) ./2;
    
    % Get the symbolic polynomial f_{i+1}
    arr_sym_fx{i} = prod(vFactors.^(vMultiplicities));
    
    % Get the degree of polynomial f_{i+1}(x)
    vDegree_fx(i) = double(feval(symengine, 'degree', (arr_sym_fx{i})));
end

% Display Polynomial f(x) in symbolic form
display(arr_sym_fx{1})

% Get the degree structure of the polynomials h_{i} where h_{i} =
% f_{i-1}(x)/f_{i}(x)
vDegree_arr_hx = diff(vDegree_fx);

% Get the degree structure of the polynomials w_{i} where w_{i} =
% h_{i-1}/h_{i}
vDegree_arr_wx = diff([vDegree_arr_hx; 0]);

% Get the multiplicity of each of the factors of f(x) by finding where the
% degree of the polynomials w_{i}(x) != 0.
vMultiplicities = find(vDegree_arr_wx ~= 0);


% Get the number of polynomials in the solution array h_{i}(x)
nPolys_arr_hx = nPolys_arr_fx - 1;


% Get the sequence of polynomials h_{i}(x) in symbolic form
sym_arr_hx = cell(nPolys_arr_hx, 1);

for i = 1 : 1 : nPolys_arr_hx
    
    sym_arr_hx{i} = arr_sym_fx{i} / arr_sym_fx{i+1};
    
end


% Initialise the arrays f_{i}(x) and h_{i}(x)
arr_fx_exact = cell(nPolys_arr_fx , 1);
arr_hx_exact = cell(nPolys_arr_hx , 1);

for i = 1:1:nPolys_arr_fx
    
    if i <= nPolys_arr_hx
        
        arr_fx_exact{i,1} = sym2poly(arr_sym_fx{i})';
        arr_hx_exact{i,1} = sym2poly(sym_arr_hx{i})';
        
    else
        
        arr_fx_exact{i,1} = 1;
        
    end
    
end

% %
% %
% %
% Add noise to the coefficients of f_{i}(x)
arr_fx_noisy = cell(nPolys_arr_fx, 1);

for i = 1 : 1 : nPolys_arr_fx
    
    arr_fx_noisy{i,1} = AddNoiseToPoly(arr_fx_exact{i},emin);
    
end


% Define an array of deconvolution methods to be used
arr_DeconvolutionMethod = {...
    'Deconvolution Separate',...
    'Deconvolution Batch',...
    'Deconvolution Batch With STLN',...
    'Deconvolution Batch Constrained',...
    'Deconvolution Batch Constrained With STLN'};


nMethods = length(arr_DeconvolutionMethod);

% Testing deconvolution
LineBreakLarge();
arr_hx = cell(nMethods,1);
arr_Error = cell(nMethods,1);

for i = 1 : 1 : nMethods
    
    % Get deconvolution method
    method_name = arr_DeconvolutionMethod{i};
    
    switch method_name
        
        case 'Deconvolution Separate'
            arr_hx{i,1} = Deconvolve_Separate(arr_fx_noisy);
            
        case 'Deconvolution Batch'
            arr_hx{i,1} = Deconvolve_Batch(arr_fx_noisy);
            
        case 'Deconvolution Batch With STLN'
            arr_hx{i,1} = Deconvolve_Batch_With_STLN(arr_fx_noisy);
            
        case 'Deconvolution Batch Constrained'
            arr_hx{i,1} = Deconvolve_Batch_Constrained(arr_fx_noisy, vMultiplicities);
            
        case 'Deconvolution Batch Constrained With STLN'
            arr_hx{i,1} = Deconvolve_Batch_Constrained_With_STLN(arr_fx_noisy, vMultiplicities);
            
        otherwise
            error('err')
            
    end
    
    fprintf([mfilename ' : ' sprintf('%s \n', method_name )]);
    
    arr_Error{i,1} = GetErrors(arr_hx{i}, arr_hx_exact);
    
end






% Plotting
nPolys_hx = size(arr_hx,1);
if (SETTINGS.PLOT_GRAPHS)
    
    figure_name = sprintf([mfilename ':' 'Deconvolution Methods Error']);
    figure('name',figure_name)
    hold on
    
    
    for i = 1:1:length(arr_DeconvolutionMethod)
        
        methodName = arr_DeconvolutionMethod{i};
        vError = arr_Error{i};
        plot(log10(vError), '-o', 'DisplayName', methodName)
        
    end
    
    
   
legend(gca,'show');
xlim([1 nPolys_hx]);
xlabel('Factor')
ylabel('log_{10} error')
hold off

end

%--------------------------------------------------------------------------
% Console writing

for i = 1:1:nMethods
    
    methodName = arr_DeconvolutionMethod{i};
    vError = arr_Error{i};
    display([mfilename ' : ' sprintf('Error %s : %e', methodName, norm(vError))]);
    
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


fullFileName = sprintf('Deconvolution/Results/Results_o_deconvolutions%s.txt',datetime('today'));
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
            bool_preproc,...
            noise,...
            method_name,...
            error_norm...
            );
        
    end

    function WriteHeader()
        
        fprintf(fileID,'DATE, EX_NUM, BOOL_PREPROC, NOISE, method_name, error_norm \n');
        
    end


end

function vErrors = GetErrors(arr_hx_comp,arr_hx_exact)
% Compare each computed h{i} with actual h_{i}
%
% % Inputs
%
% arr_hx_comp : (Array of Vectors) Each vector contains coefficients of the
% polynomial h_{i}(x)
%
% arr_hx_exact : (Array of Vectors) Each vector contains coefficients of
% the polynomial h_{i}(x)

% Get number of polynomials in the array
nPolys_hx = size(arr_hx_comp,1);

% Initialise vector to store errors
vErrors = zeros(nPolys_hx,1);

%
for i = 1:1:nPolys_hx
    
    % Get exact polynomial
    exact = arr_hx_exact{i}./ arr_hx_exact{i,1}(1,1);
    
    % Get computed polynomial
    comp = arr_hx_comp{i}./arr_hx_comp{i,1}(1,1);
    
    % Get Error
    vErrors(i) = norm(exact - comp) ./ norm(exact);
    
    
end
end
function [] = o_Deconvolution(ex_num)
% Performs a deconvolution test
%
% % Inputs
%
% ex_num : (String) Example Number
%
% % Examples
%
% o_Deconvolution('1')

% Set settings pertaining to this test

global SETTINGS
SETTINGS.PLOT_GRAPHS = true;
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-13;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 50;

% Add relevant paths
restoredefaultpath();

% Add examples folder
addpath(genpath('../Examples'));

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 

% Add that folder plus all subfolders to the path.
addpath(genpath(folder));


% Input f_{i} polynomials
x = sym('x');
y = sym('y');

% Get example
[factor_mult_arr_f] = Deconvolution_Examples_Bivariate(ex_num);

% Get vector of factors
factor = factor_mult_arr_f(:,1);

% Get vector of multiplicities of the factors
vMult = double(factor_mult_arr_f(:,2));

% Get highest power of any factor
highest_pwr = max(vMult);

% Generate polynomials f_{0}(x) ,..., f_{m}(x) = 1. Where each f_{i+1}(x) is
% the f_{i+1} = GCD(f_{i},f'_{i}).
arr_sym_fxy = cell(highest_pwr+1,1);

% Initialise a vector to store total degress of polynomials f_{i}(x,y)
vDeg_t_arr_fxy = zeros(highest_pwr+1,1);

% Initialise vectors to store degrees of f_{i}(x,y) with respect to x and y
vDeg_x_arr_fxy = zeros(highest_pwr+1,1);
vDeg_y_arr_fxy = zeros(highest_pwr+1,1);

for i = 0:1:highest_pwr
    
    % Get multiplicity of each root in f_{i+1}
    mults = ((vMult - i) + abs(vMult-i)) ./2;
    
    % Get the symbolic polynomial f_{i+1}(x,y)
    arr_sym_fxy{i+1} = prod(factor.^(mults));
    
    % Get total degree of f_{i+1}(x,y)
    vDeg_t_arr_fxy(i+1) = double(feval(symengine, 'degree', (arr_sym_fxy{i+1})));
    
    % Get the degree of f_{i+1} with respect to x
    vDeg_x_arr_fxy(i+1) = double(feval(symengine, 'degree', (arr_sym_fxy{i+1}),x));
    
    % Get the degree of f_{i+1} with respect to x
    vDeg_y_arr_fxy(i+1) = double(feval(symengine, 'degree', (arr_sym_fxy{i+1}),y));
    
end


% Get the degree structure of the polynomials h_{i}
vDegt_arr_hxy = abs(diff(vDeg_t_arr_fxy));

% Get coefficients vectors of f_{i}(x) and h_{i}(x)
nPolys_arr_fxy = size(arr_sym_fxy,1);
nPolys_arr_hxy = nPolys_arr_fxy - 1;


% Get the sequence of polynomials h_{i}(x) in symbolic form
arr_sym_hxy = cell(nPolys_arr_hxy,1);

for i = 1 : 1 : nPolys_arr_hxy
    
    arr_sym_hxy{i,1} = arr_sym_fxy{i} / arr_sym_fxy{i+1};
    
end


arr_fxy = cell(nPolys_arr_fxy, 1);
arr_hxy_exact = cell(nPolys_arr_hxy, 1);


for i = 1:1:nPolys_arr_fxy
    arr_fxy{i,1} = double(rot90(coeffs(arr_sym_fxy{i},[x,y],'All'),2));
end


for i = 1:1:nPolys_arr_hxy
    arr_hxy_exact{i,1} = double(rot90(coeffs(arr_sym_hxy{i},[x,y],'All'),2));
end

%--------------------------------------------------------------------------
% %
% %
% %

% Get number of polynomials in array f_{i}(x)
nPolys_arr_fxy = size(arr_fxy,1);

% Ensure that each f_{i}(x,y) is in a matrix of size m+1 x m+1
arr_fxy_total = cell(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy
    m = vDeg_t_arr_fxy(i);
    temp_mat = zeros(m+1,m+1);
    [nRows,nCols] = size(arr_fxy{i});
    temp_mat(1:nRows,1:nCols) = arr_fxy{i};
    arr_fxy_total{i} = temp_mat;
end

% Ensure that each h_{i}(x,y) is in a matrix of size (n+1,n+1) when
% comparing with the outputs of this deconvolution method.
arr_hxy_total = cell(nPolys_arr_hxy,1);

for i = 1:1:nPolys_arr_hxy
    n = vDegt_arr_hxy(i);
    temp_mat = zeros(n+1,n+1);
    [nRows,nCols] = size(arr_hxy_exact{i});
    temp_mat(1:nRows,1:nCols) = arr_hxy_exact{i};
    arr_hxy_total{i} = temp_mat;
end

% Deconvolution Methods

arr_method = {...
    'Separate Total',...
    'Separate Respective',...
    'Separate Both',...
    'Batch Total',...
    'Batch Respective',...
    'Batch Both'};

nMethods = length(arr_method);
arr_hxy = cell(nMethods,1);
arr_vErrors = cell(nMethods,1);
arr_error_norm = cell(nMethods,1);

for i = 1 : 1 : nMethods

    method_name = arr_method{i};
    
    switch method_name
        
        case 'Separate Total'
           
           % Intialise array
            arr_hxy_Separate_Total = cell(nPolys_arr_hxy,1);

            for j = 1 : 1 : nPolys_arr_fxy - 1

               fxy = arr_fxy{j};
               m = vDeg_t_arr_fxy(j);
               gxy = arr_fxy{j+1};
               n = vDeg_t_arr_fxy(j+1);

               arr_hxy_Separate_Total{j} = Deconvolve_Bivariate_Single_Total(fxy, gxy, m, n);
               
            end
            
            arr_hxy{i} = arr_hxy_Separate_Total;
            
        case 'Separate Respective'
            
             % Initialise array
            arr_hxy_Separate_Respective = cell(nPolys_arr_hxy, 1);

            for j = 1:1:nPolys_arr_fxy - 1

               fxy = arr_fxy{j};
               gxy = arr_fxy{j+1};
               arr_hxy_Separate_Respective{j} = Deconvolve_Bivariate_Single_Respective(fxy, gxy);

            end
            
            arr_hxy{i} = arr_hxy_Separate_Respective;
            
        case 'Separate Both'
            
            % Intialise array
            arr_hxy_Separate_Both = cell(nPolys_arr_hxy,1);

            for j = 1:1:nPolys_arr_fxy - 1

                fxy = arr_fxy{j};
                m = vDeg_t_arr_fxy(j);
                gxy = arr_fxy{j+1};
                n = vDeg_t_arr_fxy(j+1);

               arr_hxy_Separate_Both{j} = Deconvolve_Bivariate_Single_Both(fxy, gxy, m, n);
            end
            
            arr_hxy{i} = arr_hxy_Separate_Both;
        
        case 'Batch Total'
            arr_hxy{i} = Deconvolve_Bivariate_Batch_Total(arr_fxy_total, vDeg_t_arr_fxy);            
            
        case 'Batch Respective'
            
            arr_hxy{i} = Deconvolve_Bivariate_Batch_Respective(arr_fxy);
                        
        case 'Batch Both'
            arr_hxy{i} = Deconvolve_Bivariate_Batch_Both(arr_fxy, vDeg_t_arr_fxy, vDeg_x_arr_fxy, vDeg_y_arr_fxy);

    end
    
     % Get vector of error of each h_{i}(x,y) as a vector
    arr_vErrors{i} = GetError(arr_hxy{i}, arr_hxy_exact);
    arr_error_norm{i} = (norm(arr_vErrors{i}));
    
    
end


PrintToResultsFile(ex_num, arr_method, arr_error_norm);


end

function [] = PrintToResultsFile(ex_num, arr_methods, arr_error_norm)
%
% % Inputs
%
% ex_num :Example Number
%
% arr_methods : (Array of Strings)
%
% arr_error_norm : (Array of Floats)


% Get number of methods used
nMethods = length(arr_methods);

fullFileName = sprintf('Deconvolution/Results/Results_o_deconvolutions%s.txt',datetime('today'));



% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    
    for i = 1:1: nMethods
        method_name = arr_methods{i};
        error_norm = arr_error_norm{i};
        WriteNewLine(method_name, error_norm)
    end
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    for i = 1:1: nMethods
        method_name = arr_methods{i};
        error_norm = arr_error_norm{i};
        WriteNewLine(method_name, error_norm)
    end
    fclose(fileID);
    
end

    function WriteNewLine(method_name, error_norm)
        
        %
        fprintf(fileID,'%s,%s,%s,%s \n',...
            datetime('now'),...
            ex_num,...
            method_name,...
            error_norm...            
            );
        
    end

    function WriteHeader()
        
        fprintf(fileID,'DATE, EX_NUM, method_name, error_norm \n');
        
    end


end


function v_errors =  GetError(arr_hxy_comp, arr_hxy)
%
% % Inputs
%
% arr_hxy : (Array of Vectors) Vectors containing coefficients of
% polynomials h_{i}(x,y) as given exactly.
%
% arr_hxy_comp : (Array of Vectors) Vectors containing coefficients of
% polynomials h_{i}(x,y) as computed by deconvolution method.
%
% % Outputs
%
% v_errors : (Vector) Error measure for each h_{i}(x,y)


% Get number of polynomials h_{i}(x,y) in array
nPolys_arr_hxy = size(arr_hxy_comp, 1);

% Initialise vector of errors
v_errors = zeros(nPolys_arr_hxy, 1);

for i = 1:1:nPolys_arr_hxy
    
    hxy_comp_norm = arr_hxy_comp{i} ./ arr_hxy_comp{i}(1,1);
    hxy_exact = arr_hxy{i}./ arr_hxy{i}(1,1);
    
    err = abs(hxy_exact - hxy_comp_norm) ./ hxy_exact;
    
    % remove nan values
    err(isnan(err)) = [];
    
    % remove inf values
    err(isinf(err)) = [];
    
    v_errors(i) = norm(err);
   
end





end
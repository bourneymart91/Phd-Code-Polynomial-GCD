function [] = ScalingEffect_2Polys_WithCoefficients(ex_num)
% MaxMinEntriesDT(m,n_k)
%
% Each coefficient a_{i} appears in n-k+1 columns of the Sylvester matrix,
% and has two binomial coefficients in D^{-1}T(f,g).
% This experiment looks at the scaling effect of the two binomial
% coefficients of each a_{i} in each column j = 0,...,n-k.
%
% Inputs
%
% ex_num : (String) Example number
%

close all; 
clc;

% Get roots from example file
[fx, gx, ~, ~, ~] = Examples_GCD(ex_num);

% Get degree of f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Set k = 1 (Only consider first subresultant matrix.
k = 1;


% Get an array of Sylvester matrix variants
arrSylvesterMatrixType = {'DTQ', 'TQ', 'DT', 'T', 'DTQ Denominator Removed'};

% Get number of sylvester matrix variants in the array
nFormats = length(arrSylvesterMatrixType);


% For each Sylvester subresultant format
for i = 1 : 1 : nFormats
    
    
    subresultantFormat = arrSylvesterMatrixType{i};
    
    arrScaling_gx = GetScalingArray(gx, m-k, subresultantFormat);
    PlotScaling(arrScaling_gx, m-k, subresultantFormat, 'g(x)', ex_num)
    
    arrScaling_fx = GetScalingArray(fx, n-k, subresultantFormat);
    PlotScaling(arrScaling_fx, n-k, subresultantFormat, 'f(x)', ex_num)
    
    
    
    
end


end


function [arrScaling_gx] = GetScalingArray(gx, m_k, subresultantFormat)
%
% % Inputs
%
% gx : (Vector) Coefficients of polynomial g(x)
%
% m_k : (Int) number of columns in partition of k-th subresultant matrix


% Get the degree of g(x)
n = GetDegree(gx);



% For each of the coefficients b_{i}
for i = 0 : 1 : n
    
    % Get the coefficient
    bi = gx(i + 1);
    
    
    % Get vector of scaling of the coefficient a_{i} in each column of
    % the first partition of the subresultant matrix
    vScaling = bi .* GetScalingVector(n, m_k, i, subresultantFormat);
    
    % Add to array
    arrScaling_gx{i+1, 1} = abs(vScaling);
    
end

end

function [] = PlotScaling(arrScaling_gx, m_k, subresultantFormat, polyName, ex_num)
%
% arrScaling_gx : (Array of Vectors)
%
% m_k : (Int) m - k
%
% subresultantFormat : (String) 
%
% PolyName : (String) 
%
% ex_num : (String)
% 


% Get matrix where each column row i contains entries of a_{i}
gx_matrix = cell2mat(arrScaling_gx')';

% Get max in each column
vMax = max(gx_matrix);
% Get min in each column
vMin = min(gx_matrix);






% Get degree of g(x)
n = length(arrScaling_gx) - 1;

% Begin the plotting
figure_name = sprintf('%s : Scaling effect of %s in %s', mfilename, polyName, subresultantFormat );

figure('name',figure_name)
hold on

% Initialise vector x
x_vec = 0 : 1 : m_k;

for i = 0 : 1 : n
    
    % Plot values for coefficient
    
    % Set name to appear in legend
    legend_name = sprintf('Coefficient : b_{%i}',i);
    
    %plot(log10(arrScaling{i+1,1}),'DisplayName',legend_name)
    plot(x_vec, log10(arrScaling_gx{i+1,1}), 'DisplayName', legend_name)
    
end

plot(x_vec, log10(vMax),'--o', 'LineWidth',2)
plot(x_vec, log10(vMin),'--o', 'LineWidth',2)

delta = (log10(vMax) - log10(vMin));

plot(x_vec, (log10(vMin) + delta / 2), '--o', 'LineWidth',2)

xlim([0 m_k]);
%legend(gca, 'show');
ylabel('log_{10}(\Re)')
xlabel('Column index')
hold off

ResizePlot()

SaveFigure(subresultantFormat, polyName, ex_num)


end

function vScaling = GetScalingVector(m, n_k, i, sylvesterMatrixFormat)
%
% Return a vector of scaling/ coefficient multipliers for the i-th 
% coefficient of f(x), in each of the (n - k + 1) columns of S_{k}(f,g)
%
% % Inputs
%
% m : (Int) Degree of polynomial f(x)
%
% n_k : (Int) Number of columns in the k-th subresultant matrix partition
%
% i : (Int) Index of coefficient 
%
% sylvesterMatrixFormat : (String) 


vScaling = zeros(n_k + 1, 1);

switch sylvesterMatrixFormat
    case 'T'
        
        for j = 0 : 1 : n_k
            vScaling(j + 1) = nchoosek(m,i);
        end
        
    case 'DT'
        
        for j = 0 : 1 : n_k
            vScaling(j + 1) = nchoosek(m,i) ./ nchoosek(m + n_k, i + j);
        end
        
    case 'TQ'
        
        for j = 0 : 1 : n_k
            vScaling(j + 1) = nchoosek(m, i) .* nchoosek(n_k, j);
        end
        
    case 'DTQ'
        for j = 0 : 1 : n_k
            vScaling(j + 1) = nchoosek(m, i) .* nchoosek(n_k, j) ./ nchoosek(m+n_k, i + j);
        end
        
    case 'DTQ Denominator Removed'
        for j = 0 : 1 : n_k
            vScaling(j + 1) = nchoosek(m, i) .* nchoosek(n_k, j) ...
                ./ nchoosek(m+n_k, i + j) .* nchoosek(m+n_k, n_k);
        end
        
    otherwise
        error('err')
end
end


function [] = ResizePlot()

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [0.0842    0.0821    0.9018    0.9019];

end


function [] = SaveFigure(sylvester_matrix_variant, polyName, ex_num)
%
% % Inputs
%
% ex_num : (String)
%
% str : (String)
%
% sylvester_matrix_variant : (String) 

directory_name = strcat('UnivariateSylvesterFormatFigures_Coefficients/Example',(ex_num),'/');

mkdir(directory_name)

myplot = gca;

myFileName = strcat(polyName, '_', sylvester_matrix_variant);

for i = 2 : 1 : 2
    
    %saveas(h(i), [directory_name num2str(length(h) + 1 - i)], 'fig');
    saveas(myplot, [directory_name myFileName], 'fig');
    saveas(myplot, [directory_name myFileName], 'eps');
    saveas(myplot, [directory_name myFileName], 'png');
    
end


end
function [] = ScalingEffect_2Polys(m, n_k)
% ScalingEffect_2Polys(m, n_k)
%
% Each coefficient a_{i} appears in n - k + 1 columns of the k-th 
% subresultant matrix, and has two binomial coefficients in D^{-1}T(f,g).
% This experiment looks at the scaling effect of the two binomial
% coefficients of each a_{i} in each of the n - k + 1 columns.
%
% Inputs
%
% m : (Int) Degree of polynomial f(x)
%
% n_k : (Int) Number of columns in partition of k-th subresultant matrix
%
% >> ScalingEffect_2Polys(m, n_k)


close all; 
clc;

% Get an array of Sylvester matrix formats
arrSylvesterMatrixType = {'DTQ', 'TQ', 'DT', 'T', 'DTQ Denominator Removed'};

% Get number of formats in the array
nVariants = length(arrSylvesterMatrixType);

% Initialise an array to store
arrScaling = cell(m + 1, 1);

% For each subresultant matrix variant
for j = 1 : 1 : nVariants
    
    % For each of the coefficients a_{i}
    for i = 0 : 1 : m
        
        % Get format string
        subresultantFormat = arrSylvesterMatrixType{j};
        
        % Get vector of scaling of the coefficient a_{i} in each column of
        % the first partition of the subresultant matrix
        vScaling = GetScalingVector(m, n_k, i, subresultantFormat);
        
        % Add to array
        arrScaling{i+1,1} = vScaling;
        
    end
    
    
    
    % Begin the plotting
    figure_name = sprintf('%s : Scaling effect in %s',mfilename, subresultantFormat );
    figure('name',figure_name)
    hold on
    x_vec = 0 : 1 : n_k;
    
    for i = 0 : 1 : m
        % Plot values for coefficient
        
        % Set name to appear in legend
        legend_name = sprintf('Coefficient : a_{%i}',i);
        
        %plot(log10(arrScaling{i+1,1}),'DisplayName',legend_name)
        plot( x_vec, (arrScaling{i+1,1}), 'DisplayName', legend_name)
    end
    
    xlim([0 n_k]);
    legend(gca, 'show');
    ylabel('log_{10}(\Re)')
    xlabel('Column index')
    hold off
    
    
    
    
    
    
    
    
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [0.0842    0.0821    0.9018    0.9019];
    
    % Set line width
    
    
    
    
    
    %ax.Position = [left bottom ax_width ax_height];
    
    grid on
    box on
    
    hLines = findobj(gca,'Type','line');

    for k = 1:length(hLines)
        activeLine = hLines(k);
        
        set(activeLine,'LineWidth',6)
        set(activeLine,'markers',1)
        set(activeLine,'Marker','square')
        set(activeLine, 'MarkerEdgeColor','black')
        set(activeLine,'MarkerFaceColor','black')
        
    end
    
end



end



function vScaling = GetScalingVector(m, n_k, i, sylvesterMatrixFormat)
%
% % Inputs
%
% m : (Int)
%
% n_k : (Int)
%
% i : (Int)
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
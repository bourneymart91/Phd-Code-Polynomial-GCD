function [] = SylvesterMatrixHeatMap_2Polys(m, n, k)
%
% % Inputs
%
% m : (Int) Degree of polynomial f(x)
%
% n : (Int) Degree of polynomial g(x)
%
% k : (Int) Index of k-th subresultant matrix
%
% % Outputs
%
% % Example
%
% >> SylvesterMatrixHeatMap(5, 15, 4)

close all; 
clc;


% Initialise two polynomials with coefficients all set to '1'.
fx = ones(m + 1, 1);
gx = ones(n + 1, 1);

% Get an array of subresultant matrix variants
arrSubresultantFomat = {'T','DT','TQ','DTQ', 'DTQ Denominator Removed'};

% Get number of variants considered
nFormats = length(arrSubresultantFomat);

% For each variant of subresultant matrix
for i = 1 : 1 : nFormats
    
    
    subresultantVariant = arrSubresultantFomat{i};
    
    
    switch subresultantVariant
        
        case 'T'
            
            T1 = BuildT1(fx, n - k);
            T2 = BuildT1(gx, m - k);
            
            Sk = [T1 T2];
            
        case 'DT'
            
            D = BuildD_2Polys(m, n - k);
            T1 = BuildT1(fx, n - k);
            T2 = BuildT1(gx, m - k);
            
            Sk = D*[T1 T2];
            
        case 'TQ'
            
            T1 = BuildT1(fx, n - k);
            T2 = BuildT1(gx, m - k);
            Q = BuildQ_2Polys(m, n, k);
            
            Sk = [T1 T2] * Q;
            
        case 'DTQ'
            
            DT1Q1 = BuildDT1Q1(fx, n - k);
            DT2Q2 = BuildDT1Q1(gx, m - k);
            
            Sk = [DT1Q1 DT2Q2];
            
        case 'DTQ Denominator Removed'
            
            
            DT1Q1 = BuildDT1Q1_Rearranged_RemovedDenom(fx, n - k);
            DT2Q2 = BuildDT1Q1_Rearranged_RemovedDenom(gx, m - k);
            
            Sk = [DT1Q1 DT2Q2];
            
        otherwise
            error('%s is not a valid format', subresultantVariant)
            
    end
    
    figure_name = sprintf('%s : Heat Map', subresultantVariant);
    figure('Name',figure_name)
    
    Sk_rounded = log10(Sk);
    
    hm = heatmap((Sk_rounded));
    
    
end








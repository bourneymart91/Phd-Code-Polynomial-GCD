function [] = ConditionNumbersExperiments_2Polys(ex_num)
%
% % Inputs
%
% ex_num : (String) Example number


close all; 
clc;

% Get coefficients of f(x) and g(x)
[fx, gx] = Examples_GCD_FromCoefficients(ex_num);

% Get degree of polynomial f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Get an array of formats of subresultant matrix
arrFormat = {'T','DT','TQ','DTQ', 'DTQ Denominator Removed'};

% Get number of formats
nFormats = length(arrFormat);

% Initialise arrays
arrConditionVec = cell(nFormats, 1);
arrConditionVec_P1 = cell(nFormats, 1);
arrConditionVec_P2 = cell(nFormats, 1);

for i = 1:1:nFormats
    
    strFormat = arrFormat{i};
    
    vConditionSk = zeros(min(m,n), 1);
    vCondition_P1 = zeros(min(m,n), 1);
    vCondition_P2 = zeros(min(m,n), 1);
    
    for k = 1:1:min(m,n)
        
        
        
        Sk = BuildSubresultant(fx,gx,k, strFormat);
        
        vConditionSk(k) = cond(Sk);
        
        vCondition_P1(k) = cond(Sk(:, 1 : n - k + 1));
        vCondition_P2(k) = cond(Sk(:, n - k + 2 : end));
        
    end
    
    arrConditionVec_P1{i} = vCondition_P1;
    arrConditionVec_P2{i} = vCondition_P2;
    arrConditionVec{i} = vConditionSk;
    
end


x_vec = 1:1:min(m,n);
x_low = 1;
x_high = min(m,n);


% Plot Condition of Both Partitions
figure('Name','Condition Sk')
hold on
for i = 1:1:nFormats
    
    vConditionSk = arrConditionVec{i};
    FormatName = arrFormat{i};
    plot(x_vec, log10(vConditionSk),'DisplayName',FormatName,'LineWidth',2);
    xlim([x_low, x_high])
    
end
legend(gca,'show');
hold off


% Plot Condition of First Partition

figure('Name', 'Condition Cf')
hold on
for i = 1:1:nFormats
    
    vConditionSk = arrConditionVec_P1{i};
    FormatName = arrFormat{i};
    plot(x_vec, log10(vConditionSk),'DisplayName',FormatName,'LineWidth',2);
    xlim([x_low, x_high])
    
end
legend(gca,'show');
hold off


% Plot Condition of Second Partition

figure('Name','Condition Cg')
hold on
for i = 1:1:nFormats
    
    vConditionSk = arrConditionVec_P2{i};
    FormatName = arrFormat{i};
    plot(x_vec, log10(vConditionSk),'DisplayName',FormatName,'LineWidth',2);
    xlim([x_low, x_high])
end
legend(gca,'show');
hold off
end


function [Sk] = BuildSubresultant(fx, gx, k, strFormat)
% BuildSubresultant_2Polys(fx,gx,k)
%
% This function builds the k-th Sylvester subresultant matrix S_{k}(f,g),
% in the Bernstein Basis.
%
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x) in Bernstein Basis
%
% gx : (Vector) Coefficients of polynomial g(x) in Bernstein Basis
%
% k:  (Int) Index of subresultant S_{k} to be built
%
% % Outputs
%
% Sk : (Matrix) The kth Sylvester subresultant matrix S_{k}(f,g)


% Get the degree of the polynomial f(x)
m = GetDegree(fx);

% Get the degree of the polynomial g.
n = GetDegree(gx);


switch strFormat
    
    case 'T'
        
        T = BuildT_2Polys(fx, gx, k);
        Sk = T;
        
    case 'DT'
        
        D = BuildD_2Polys(m, n - k);
        T = BuildT_2Polys(fx, gx, k);
        Sk = D*T;
        
    case 'DTQ'
        
        D = BuildD_2Polys(m, n - k);
        T = BuildT_2Polys(fx, gx, k);
        Q = BuildQ_2Polys(m, n, k);
        Sk = D*T*Q;
        
    case 'TQ'
        
        T = BuildT_2Polys(fx, gx, k);
        Q = BuildQ_2Polys(m, n, k);
        
        Sk = T*Q;
        
    case 'DTQ Denominator Removed'
        
        DT1Q1 = BuildDT1Q1_Rearranged_RemovedDenom(fx, n - k);
        DT2Q2 = BuildDT1Q1_Rearranged_RemovedDenom(gx, m - k);
        
        Sk = [DT1Q1 DT2Q2];
        
        
        
    case 'DTQ Rearranged' % Build based on generating individual elements
        
        % Build First Partition
        D1T1Q = BuildDT1Q1(fx, n - k);
        
        % Build Second Partition
        D2T2Q = BuildDT1Q1(gx, m - k);
        
        % Build Sylvester Matrix by concatenation of matrices A and B.
        Sk = [D1T1Q D2T2Q];
        
    otherwise
        error('SETTINGS.SYLVESTER_MATRIX_VARIANT must be one of the following *T or *DT or *DTQ or *TQ or *DTQ Denominator Removed or *DTQ Rearranged')
end

end


function [] = CoefficientHeatMap_3Polys(m , n, o, k)
% 
%
% % Input
%
% m : (Int) Degrree of f(x,y)
%
% n : (Int) Degree of g(x,y)
%
% o : (Int) Degree of h(x,y)
%
% k : (Int)
%
% >> CoefficientHeatMap_3Polys(10, 5, 3, 2)

close all; clc;

fxy = ones(m + 1, m + 1);
gxy = ones(n + 1, n + 1);
hxy = ones(o + 1, o + 1);


arrSylvesterFormats = {...
    'T', ... 
    'DT', ...
    'TQ', ...
    'DTQ' ...
    };

nFormats = length(arrSylvesterFormats);

arrSk = cell(nFormats,1);

for i = 1:1:nFormats
   
    arrSk{i} = BuildSylvesterMatrix_3Polys(fxy, gxy, hxy, m, n, o, k, arrSylvesterFormats{i});
    
end


for i = 1:1:nFormats
   
    format_name = arrSylvesterFormats{i};
    Sk = arrSk{i};
    
    PlotHeatMap(Sk, format_name);
    
end

end




function [] = PlotHeatMap(Sk, method_name)
%
% % Inputs
%
% Sk : (Matrix)
%
% method_name : (String)

Sk = abs(Sk);
Sk = log10(Sk);


figure('name',method_name)

h = heatmap(Sk);

h.GridVisible = 'off';



end


function Sk = BuildSylvesterMatrix_3Polys(fxy, gxy, hxy, m, n, o, k, sylvester_matrix_variant)
% BuildSylvesterMatrix(fxy,gxy,m,n,k)
%
% Build the Sylvester subresultant matrix S_{k}(f,g)
% 
%
% Inputs
%
% fxy : (Matrix)
%
% gxy : (Matrix)
%
% hxy : (Matrix) Coefficients of the polynomial f(x,y), g(x,y) and
% h(x,y)
%
% m : (Int) Degree of f(x,y)
%
% n : (Int) Degree of g(x,y)
%
% o : (Int) degree of h(x,y)
%
% k : Index of the Sylvester subresultant matrix S_{k} to be constructed.
%
% sylvester_matrix_variant : (String)


switch sylvester_matrix_variant
    case 'T'
        
        
         % Build the matrix T
        T1_f = BuildT1(fxy, m, n - k);
        
        T2_f = BuildT1(fxy, m, o - k);
        
        T3_g = BuildT1(gxy, n, m - k);
        
        T4_h = BuildT1(hxy, o, m - k);
        
        diagonal = blkdiag(T1_f, T2_f);
        column = [T3_g; T4_h];
        Sk = [diagonal column];
        
    case 'DT'
        
        
        % Build the matrix D
        D = BuildD_3Polys(m, n, o, k);
        %D1 = BuildD_2Polys(m,n-k);
        %D2 = BuildD_2Polys(m,o-k);
        
        % Build the matrix T
        T1_f = BuildT1(fxy, m, n - k);
        T2_f = BuildT1(fxy, m, o - k);
        T3_g = BuildT1(gxy, n, m - k);
        T4_h = BuildT1(hxy, o, m - k);
        
        diagonal = blkdiag(T1_f, T2_f);
        column = [T3_g; T4_h];
        T = [diagonal column];

        % Build the matrix DT
        Sk = D*T;
        
    case 'DTQ'
        
         % Build the matrix D
        D = BuildD_3Polys(m, n, o, k);
        %D1 = BuildD_2Polys(m,n-k);
        %D2 = BuildD_2Polys(m,o-k);
        
        % Build the matrix T
        T1_f = BuildT1(fxy, m, n - k);
        T2_f = BuildT1(fxy, m, o - k);
        T3_g = BuildT1(gxy, n, m - k);
        T4_h = BuildT1(hxy, o, m - k);
        
        diagonal = blkdiag(T1_f, T2_f);
        column = [T3_g; T4_h];
        T = [diagonal column];
        
        Q = BuildQ_3Polys(m, n, o, k);
        
        
        Sk = D * T * Q;
    
    case 'TQ'
        
        % Build the matrix T
        T1_f = BuildT1(fxy, m, n - k);
        T2_f = BuildT1(fxy, m, o - k);
        T3_g = BuildT1(gxy, n, m - k);
        T4_h = BuildT1(hxy, o, m - k);
        
        diagonal = blkdiag(T1_f, T2_f);
        column = [T3_g; T4_h];
        T = [diagonal column];
        
        Q = BuildQ_3Polys(m, n, o, k);
        
        
        Sk = T * Q;
        
    case 'DTQ Denominator Removed'
        
        error('Not Yet Complete')
        
    otherwise
        
        error('err')
        
end


end

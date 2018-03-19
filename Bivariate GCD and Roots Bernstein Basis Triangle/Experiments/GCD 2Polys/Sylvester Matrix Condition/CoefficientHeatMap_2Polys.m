function [] = CoefficientHeatMap_2Polys(m , n, k)
% Plot heatmap of scaling effect of 
%
% % Input
%
% m : (Int) Degree of polynomial f(x,y)
%
% n : (Int) Degree of polynomial g(x,y)
%
% k : (Int) Index of subresultant matrix
%
% >> CoefficientHeatMap(10, 5, 3)

close all; 
clc;

% 
fxy = ones(m + 1, m + 1);
gxy = ones(n + 1, n + 1);



arrSylvesterFormats = {...
    'T', ... 
    'DT', ...
    'TQ', ...
    'DTQ', ...
    'DTQ Denominator Removed' ...
    };
% Get number of sylvester subresultant matrix variants
nFormats = length(arrSylvesterFormats);

arrSk = cell(nFormats,1);

for i = 1 : 1 : nFormats
   
    arrSk{i} = BuildSylvesterMatrix_2Polys(fxy, gxy, m, n, k, arrSylvesterFormats{i});
    
end


for i = 1 : 1 : nFormats
   
    format_name = arrSylvesterFormats{i};
    Sk = arrSk{i};
    
    PlotHeatMap(Sk, format_name);
    
end

end


function Sk = BuildSylvesterMatrix_2Polys(fxy, gxy, m, n, k, sylvester_format)
% BuildSylvesterMatrix(fxy,gxy,m,n,k)
%
% Build the Sylvester subresultant matrix S_{k}(f,g)
% 
%
% Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% m : (Int) Total degree of polynomial f(x,y)
%
% n : (Int) Total degree of polynomial g(x,y)
%
% k : (Int) Index of the Sylvester subresultant matrix S_{k} to be constructed.
%
% % Outputs
%
% Sk : (Matrix) kth Sylvester subresultant matrix S_{k}(f,g)




switch sylvester_format
    case 'T'

        % Build the matrix T_{n-k}(f(x,y))
        T1_fx = BuildT1(fxy, m, n - k);

        % Build the matrix T_{m-k}(g(x,y))
        T1_gx = BuildT1(gxy, n, m - k);

        Sk = [T1_fx T1_gx] ;
        
    case 'DT'
        % Build the diagonal matrix D^{-1}
        D = BuildD_2Polys(m, n - k);

        % Build the matrix T_{n-k}(f(x,y))
        T1_fx = BuildT1(fxy, m, n - k);

        % Build the matrix T_{m-k}(g(x,y))
        T1_gx = BuildT1(gxy, n, m - k);

       
        Sk = D*[T1_fx T1_gx];
        
    case 'DTQ'
        % Build the diagonal matrix D^{-1}
        D = BuildD_2Polys(m, n - k);

        % Build the matrix T_{n-k}(f(x,y))
        T1_fx = BuildT1(fxy, m, n - k);

        % Build the matrix T_{m-k}(g(x,y))
        T1_gx = BuildT1(gxy, n, m - k);

        Q = BuildQ_2Polys(m, n, k);
        
        Sk = D*[T1_fx T1_gx]*Q;
    
    case 'TQ'

        % Build the matrix T_{n-k}(f(x,y))
        T1_fx = BuildT1(fxy, m, n - k);

        % Build the matrix T_{m-k}(g(x,y))
        T1_gx = BuildT1(gxy, n, m - k);

        Q = BuildQ_2Polys(m, n, k);
        
        Sk = [T1_fx T1_gx]*Q;
        
    case 'DTQ Denominator Removed'
        
        % Build the diagonal matrix D^{-1}
        D = BuildD_2Polys(m, n - k);

        % Build the matrix T_{n-k}(f(x,y))
        T1_fx = BuildT1(fxy, m, n - k) .* nchoosek(m + n - k, n - k);

        % Build the matrix T_{m-k}(g(x,y))
        T1_gx = BuildT1(gxy, n, m - k) .* nchoosek(m + n - k , m - k);

        Q = BuildQ_2Polys(m, n, k);
        
        Sk = D*[T1_fx T1_gx]*Q;
        
    otherwise
        error('err')
        
end


end


function [] = PlotHeatMap(Sk, method_name)
%
% % Inputs
%
% Sk : (Matrix)
%
% method_name : (String)


figure('name',method_name)


% Note - neccesary to take absolute values since one partition is
% negative
Sk_rounded = log10(abs(Sk));

hm = heatmap((Sk_rounded));
    
    
    



end
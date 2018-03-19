function [] = SylvesterMatrixHeatMap_3Polys_Changing_k(m, n, o)
%
% % Inputs
%
% m : (Int) Degree of f(x)
%
% n : (Int) Degree of g(x)
%
% o : (Int) Degree of h(x)
%
% k : (Int) index of Sylvester subresultant matrix
%
% % Outputs
%
%
% % Example
%
% >> SylvesterMatrixHeatMap_3Polys(5, 15, 7)

close all; clc;



fx = ones(m + 1, 1);
gx = ones(n + 1, 1);
hx = ones(o + 1, 1);

arrSubresultantFomat = {'T'};

nFormats = length(arrSubresultantFomat);

for k = 1:1:min([m,n,o])
    
    for i = 1:1:nFormats
        
        subresultantFomat = arrSubresultantFomat{i};
        
        
        switch subresultantFomat
            
            case 'T'
                
                Tf = BuildT1(fx, n - k);
                Tf2 = BuildT1(fx, o - k);
                Tg = BuildT1(gx, m - k);
                Th = BuildT1(hx, m - k);
                
                Sk = [blkdiag(Tf, Tf2) [Tg ; Th]];
                
            case 'DT'
                
                
                D = BuildD_3Polys(m, n - k, o - k);
                
                Tf = BuildT1(fx, n - k);
                Tf2 = BuildT1(fx, o - k);
                Tg = BuildT1(gx, m - k);
                Th = BuildT1(hx, m - k);
                
                Sk = D * [blkdiag(Tf, Tf2) [Tg ; Th]];
                
            case 'TQ'
                
                Tf = BuildT1(fx, n - k);
                Tf2 = BuildT1(fx, o - k);
                Tg = BuildT1(gx, m - k);
                Th = BuildT1(hx, m - k);
                Q = BuildQ_3Polys(m, n, o, k);
                
                Sk = [blkdiag(Tf, Tf2) [Tg ; Th]] * Q;
                
                
            case 'DTQ'
                
                
                D = BuildD_3Polys(m, n - k, o - k);
                
                Tf = BuildT1(fx, n - k);
                Tf2 = BuildT1(fx, o - k);
                Tg = BuildT1(gx, m - k);
                Th = BuildT1(hx, m - k);
                
                Q = BuildQ_3Polys(m, n, o, k);
                
                Sk = D * [blkdiag(Tf, Tf2) [Tg ; Th]] * Q;
                
                
            case 'DTQ Denominator Removed'
                
                error('not completed')
                
                
            otherwise
                error('%s is not a valid format', subresultantFomat)
                
        end
        
        figure_name = sprintf('%s : Heat Map', subresultantFomat);
        figure('Name',figure_name)
        
        Sk_rounded = log10(Sk);
        
        hm = heatmap((Sk_rounded));
        colormap(hm, 'hot')
        
    end
end







function [] = SylvesterMatrixHeatMap_3Polys(m, n, o, k)
% This function plots the heatmap of the combinatorial terms in each 
% variant of the k-th subresultant matrix. 
%
% % Inputs
%
%
% m : (Int) The degree of the polynomial f(x)
%
% n : (Int) The degree of the polynomial g(x)
%
% o : (Int) The degree of the polynomial h(x)
%
% k : (Int) index of kth Sylvester subresultant matrix
%
%
% % Outputs
%
%
% This function generates the heatmaps of the combinatorial terms in each 
% of the variants of the k-th subresultant matrix.  
%
%
% % Example
%
%
% >> SylvesterMatrixHeatMap_3Polys(5, 15, 7, 3, true)




% Initialise vectors of coefficients of f(x), g(x) and h(x)
fx = ones(m + 1, 1);
gx = ones(n + 1, 1);
hx = ones(o + 1, 1);


% Initialise an array of the variants of the k-th subresultant matrix
arrSubresultantMatrixVariant = {'T', 'DT', 'TQ', 'DTQ'};

% Get number of variants
nVariants = length(arrSubresultantMatrixVariant);

% Only consider the 2 x 3 partitioned subresultant matrix
SYLVESTER_EQUATIONS = '2';

% For each variant of subresultant matrix
for i = 1 : 1 : nVariants
    
    subresultantVariant = arrSubresultantMatrixVariant{i};
    
    
    switch subresultantVariant
        
        case 'T'
            
            switch SYLVESTER_EQUATIONS
                
                case '2'
                    T = BuildT_3Polys_2Eqns(fx, gx, hx, k);
                    
                case '3'
                    T = BuildT_3Polys_3Eqns(fx, gx, hx, k);
                    
                otherwise
                    error('err')
                    
            end
            
            
            Sk = T;
            
        case 'DT'
            
            switch SYLVESTER_EQUATIONS
                
                case '2'
                    
                    D = BuildD_3Polys_2Eqns(m, n - k, o - k);
                    T = BuildT_3Polys_2Eqns(fx, gx, hx, k);
                    
                case '3'
                    D = BuildD_3Polys_3Eqns(m, n, o, k);
                    T = BuildT_3Polys_3Eqns(fx, gx, hx, k);
                otherwise
                    error('err')
            end
            
            Sk = D*T;
            
        case 'TQ'
            
            switch     SYLVESTER_EQUATIONS
                case '2'
                    T = BuildT_3Polys_2Eqns(fx, gx, hx, k);
                    Q = BuildQ_3Polys(m, n, o, k);
                case '3'
                    T = BuildT_3Polys_3Eqns(fx, gx, hx, k);
                    Q = BuildQ_3Polys(m,n,o,k);
                otherwise
                    error('err')
                    
                    
            end
            
            Sk = T*Q;
            
            
        case 'DTQ'
            
            switch SYLVESTER_EQUATIONS
                case '2'
                    D = BuildD_3Polys_2Eqns(m, n-k, o-k);
                    T = BuildT_3Polys_2Eqns(fx, gx, hx, k);
                    Q = BuildQ_3Polys(m, n, o, k);
                    
                case '3'
                    D = BuildD_3Polys_3Eqns(m,n,o,k);
                    T = BuildT_3Polys_3Eqns(fx, gx, hx, k);
                    Q = BuildQ_3Polys(m, n, o, k);
                    
                otherwise
                    error('err')
            end
            
            Sk = D*T*Q;
            
            
            
        case 'DTQ Denominator Removed'
            
            error('not completed')
            
            
        otherwise
            error('%s is not a valid format', subresultantVariant)
            
    end
    
    figure_name = sprintf('%s : Heat Map', subresultantVariant);
    figure('Name',figure_name)
    
    % Note - neccesary to take absolute values since one partition is
    % negative
    Sk_rounded = log10(abs(Sk));
    
    hm = heatmap((Sk_rounded));
    
end


end




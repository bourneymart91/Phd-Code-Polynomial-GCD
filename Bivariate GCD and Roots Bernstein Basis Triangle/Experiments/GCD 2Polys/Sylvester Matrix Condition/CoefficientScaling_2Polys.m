% This experiment looks at the scaling of the coefficients of the
% polynomial f(x) in the Sylvester matrix S_{k}(f,g)

function [] = CoefficientScaling(m, n, k)
%
%
% % Inputs
%
% m : (Int) Degree of f(x,y)
%
% n : (Int) Degree of g(x,y)
%
% k : (Int) Degree of d(x,y)
%
% >> CoefficientScaling(10, 5, 3)


% Produce array of Sylvester subresultant formats
arrSylvesterFormats = {'T', 'DT', 'TQ', 'DTQ', 'DTQ Denominator Removed'};
nSylvesterFormats = length(arrSylvesterFormats);

% Get scaling of each coefficient
arrT = cell( m + 1, m + 1);
arrDT = cell( m + 1, m + 1);
arrTQ = cell( m + 1, m + 1);
arrDTQ = cell( m + 1, m + 1);
arrDTQ_DenominatorRemoved = cell( m + 1, m + 1);


for k1 = 0 : 1  : m
    
    for i1 = k1 : -1 : 0
        
        i2 = k1 - i1;
        
        % For each coefficient a_{i_{1}, i_{2}}
        
        for i = 1:1:nSylvesterFormats
            arrT{i1 + 1, i2 + 1} = GetScaling(m, n, k, i1, i2, 'T');
            arrDT{i1 + 1, i2 + 1} = GetScaling(m, n, k, i1, i2, 'DT');
            arrTQ{i1 + 1, i2 + 1} = GetScaling(m, n, k, i1, i2, 'TQ');
            arrDTQ{i1 + 1, i2 + 1} = GetScaling(m, n, k, i1, i2, 'DTQ');
            arrDTQ_DenominatorRemoved{i1 + 1, i2 + 1} = GetScaling(m, n, k, i1, i2, 'DTQ Denominator Removed');
        end
        
    end
    
    
end

% Plot scaling of each of the coefficients
PlotValues(arrT, m, n, k, 'T');
PlotValues(arrDT, m, n, k, 'DT');
PlotValues(arrTQ, m , n, k, 'TQ');
PlotValues(arrDTQ, m, n, k, 'DTQ');
PlotValues(arrDTQ_DenominatorRemoved, m, n, k, 'DTQ Denominator Removed');


end




function [] = PlotValues(arrT, m, n, k, sylvester_format)
%
% % Inputs
%
% arrT : (Array of matrices) Each matrix contains the scaling of coefficient
%   a_{i1,i2} in each of the j_{1},j_{2} columns
%
% m : (Int) Degree of f(x,y)
%
% N : (Int) Degree of g(x,y)
%
% k : (Int) Degree of d(x,y)
%
% sylvester_format : (string) 


figure('name',sylvester_format)

nNonZeros_vxy = nchoosek(n - k + 2, 2);





hold on
for k1 = 0 : 1  : m
    
    for i1 = k1 : -1 : 0
        
        i2 = k1 - i1;
        
        temp_vec = GetAsVector(arrT{i1 + 1, i2 + 1});
        temp_vec = temp_vec(1: nNonZeros_vxy);
        temp_vec = log10(temp_vec);
        
        x_vec = 1:1:length(temp_vec);
        plot_name = sprintf('a_{%i,%i}',i1,i2);
        plot(x_vec, temp_vec,'-s', 'DisplayName', plot_name);
        
    end
end
legend(gca,'show');

hold off


end


function [scaling_matrix] = GetScaling(m, n, k, i1, i2, method_name)
%
% % Inputs
%
% m : (Int) Degree of f(x,y)
%
% n : (Int) Degree of g(x,y)
%
% k : (Int) Degree of d(x,y)
%
% i1 : (Int)
%
% i2 : (Int) 
%
% method_name : (String)
%
% % Outputs
%
% scaling_matrix : (Matrix)

scaling_matrix = zeros(n - k + 1, n - k + 1);


for k2 = 0 : 1 : n - k
    
    for j1 = k2: -1 : 0
        
        j2 = k2 - j1;
        
        % For each column
        
        
        switch method_name
            
            case 'T'
                scaling_matrix(j1 + 1, j2 + 1) = Trinomial(m, i1, i2);
                
            case 'DT'
                scaling_matrix(j1 + 1, j2 + 1) = ...
                    Trinomial(m, i1, i2)...
                    ./ ...
                    Trinomial(m + n - k, i1 + j1, i2 + j2);
                
            case 'TQ'
                scaling_matrix(j1 + 1, j2 + 1) = ...
                    Trinomial(m, i1, i2) ...
                    * ...
                    Trinomial(n - k, j1, j2);
                
            case 'DTQ'
                scaling_matrix(j1 + 1, j2 + 1) = ...
                    Trinomial(m, i1, i2) ...
                    * ...
                    Trinomial(n - k, j1, j2) ...
                    ./ ...
                    Trinomial(m + n - k, i1 + j1, i2 + j2);
                
            case 'DTQ Denominator Removed'
                
                scaling_matrix(j1 + 1, j2 + 1) = ...
                    nchoosek(i1 + j1, i1) ...
                    * ...
                    nchoosek(i2 + j2, i2) ...
                    / ...
                    nchoosek(m + n - k - i1 - i2 - j1 - j2, m - i1 - i2) ...
                    .* ...
                    nchoosek(m + n - k, n - k);
                
                
        end
        
    end
end

end

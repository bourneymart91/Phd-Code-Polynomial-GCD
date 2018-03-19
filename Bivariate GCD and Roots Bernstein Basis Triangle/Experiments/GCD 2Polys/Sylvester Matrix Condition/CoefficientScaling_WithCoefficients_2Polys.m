% This experiment looks at the scaling of the coefficients of the
% polynomial f(x) in the Sylvester matrix S_{k}(f,g)

function [] = CoefficientScaling_WithCoefficients_2Polys(ex_num)
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
% >> CoefficientScaling_WithCoefficients_2Polys('1')



% %
% %
% Get two polynomials f(x,y) and g(x,y) from example file.
[fxy, gxy, uxy_exact, vxy_exact, dxy_exact, m, n, t_exact] = Examples_GCD(ex_num);

m = GetDegree_Bivariate(fxy);
n = GetDegree_Bivariate(gxy);
k = 1;

% Produce array of Sylvester subresultant formats
arrSylvesterFormats = {'T', 'DT', 'TQ', 'DTQ', 'DTQ Denominator Removed'};
nSylvesterFormats = length(arrSylvesterFormats);


for i = 1:1:nSylvesterFormats
     
    sylvesterFormat = arrSylvesterFormats{i};
    GetEntries(fxy, n, k, sylvesterFormat, 'f(x)');
    GetEntries(gxy, m, k, sylvesterFormat, 'g(x)');
    
end







end

function [] = GetEntries(fxy, n, k, sylvester_format, polyName)

m = GetDegree_Bivariate(fxy);

for k1 = 0 : 1  : m
    
    for i1 = k1 : -1 : 0
        
        i2 = k1 - i1;
        
        ai1i2 = fxy(i1+1,i2+1);
        
        % For each coefficient a_{i_{1}, i_{2}}
        arrT{i1 + 1, i2 + 1} = ai1i2 .* GetScaling(m, n, k, i1, i2, sylvester_format);
        
    end
    
    
end

% Plot scaling of each of the coefficients
PlotValues(arrT, m, n, k, sylvester_format, polyName);


end


function [] = PlotValues(arrT, m, n, k, sylvester_format, polyName)
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


figure('name',strcat(sylvester_format,'_', polyName));

nColumns = nchoosek(n - k + 2, 2);





hold on
for k1 = 0 : 1  : m
    
    for i1 = k1 : -1 : 0
        
        i2 = k1 - i1;
        
        temp_vec = GetAsVector(arrT{i1 + 1, i2 + 1});
        temp_vec = temp_vec(1: nColumns);
        temp_vec = log10(abs(temp_vec));
        
        x_vec = 1:1:length(temp_vec);
        plot_name = sprintf('a_{%i,%i}',i1,i2);
        plot(x_vec, temp_vec,'-s', 'DisplayName', plot_name);
        
    end
end

%legend(gca,'show');

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

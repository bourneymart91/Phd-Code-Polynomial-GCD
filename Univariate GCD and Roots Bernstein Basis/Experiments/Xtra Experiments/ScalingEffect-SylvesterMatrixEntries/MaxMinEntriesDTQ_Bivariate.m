function [] = MaxMinEntriesDTQ_Bivariate(m1,m2,n1_k1,n2_k2)
% MaxMinEntriesDTQ_Bivariate(m, n_k)
%
% Each coefficient a_{i} appears in (n1-k1+1) \times (n2 - k2 + 1) columns 
% of the Sylvester matrix, and has three binomial coefficients in 
% D^{-1}T(f,g)Q.
% This experiment looks at the scaling effect of the three binomial
% coefficients of each a_{i} in each column j = 0,...,n-k.
%
% Inputs
%
% m1 : (Int) Degree of polynomial f(x,y) with respect to x
%
% m2 : (Int) Degree of polynomail f(x,y) with respect to y
%
% n1_k1 : (Int)
%
% n2_k2 : (Int)
%
% Example
%
% >> MaxMinEntriesDTQ_Bivariate(7, 5, 2, 2)



figure_name = sprintf('%s : Scaling effect',mfilename);
figure('name',figure_name)
hold on

% Initialise a matrix to store entries of the k1k2-th subresultant matrix
a = zeros(m1 + 1, m2 + 1, n1_k1 + 1, n2_k2 + 1);

% For each coefficient a_{i} of polynomial f(x)
for i1 = 0:1:m1
    for i2 = 0:1:m2
        
        % for each column (j_{1}, j_{2})
        for j1 = 0:1:n1_k1
            for j2 = 0:1:n2_k2
                
                
                a(i1+1,i2+1,j1+1 ,j2 +1) = ...
                    nchoosek(m1,i1) ...
                    * nchoosek(m2,i2) ...
                    * nchoosek(n1_k1,j1) ...
                    * nchoosek(n2_k2,j2) ...
                    ./ nchoosek(m1+n1_k1,i1+j1) ...
                    ./ nchoosek(m2+n2_k2,i2+j2);
                
            end
        end      
    end
end

% Plotting
legend(gca,'show');
ylabel('Scaling of coefficient')
xlabel('Column index')
hold off

figure()
hold on
for i1 = 0:1:m1
    for i2 = 0:1:m2
        
        surface = squeeze(a(i1+1, i2+1, :, :));
        surf(log10(surface))
        
    end
end
hold off

end
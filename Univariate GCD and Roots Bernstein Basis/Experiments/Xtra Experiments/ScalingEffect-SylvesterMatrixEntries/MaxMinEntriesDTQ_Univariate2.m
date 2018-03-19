function [] = MaxMinEntriesDTQ_Univariate2(m, max_n_k)
% MaxMinEntriesDTQ_Univariate(m,n_k)
%
% Each coefficient a_{i} appears in n-k+1 columns of the Sylvester matrix,
% and has three binomial coefficients in D^{-1}T(f,g)Q.
% This experiment looks at the scaling effect of the three binomial
% coefficients of each a_{i} in each column j = 0,...,n-k.
%
% Inputs
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial v(x)
%
% >> MaxMinEntriesDTQ_Univariate(m,n_k)




vScalar = zeros(m+1, max_n_k+1, max_n_k + 1);
vProduct = zeros(1,m+1);
vSum = zeros(1,m+1);

% for each a_{i}
for i = 0:1:1
    for n_k = 1:1:max_n_k
    for j = 0:1:n_k
        vScalar(i+1,n_k,j+1) = (nchoosek(i+j,i) ...
            * nchoosek(m+n_k-(i+j),m-i)) ...
            / nchoosek(m+n_k,m);
        
    end
    end
    
    
end

% % Plot
figure_name = sprintf('%s : Scaling effect',mfilename);
figure('name',figure_name)
hold on
for i = 0:1:m
    name = sprintf('Coefficient : a_{%i}',i);
    %vProduct(i+1) = prod(vScalar(i+1,:));
    %vSum(i+1) = sum(vScalar(i+1,:));
    surf(log10(squeeze((vScalar(i+1,:,:)))),'DisplayName',name,'FaceColor',rand(1,3));
end
legend(gca,'show');
zlabel('Scaling of coefficient')
xlabel('Column index')
ylabel('n-k')
hold off


% % Take all entries of C(f) containing each a_{i}, and get the product of
% % all entries containining a_{i}. \prod_{j=0,...,n-k}
% figure_name = sprintf('%s : Product of coefficient occurences in C(f)',mfilename);
% figure('name',figure_name)
% hold on
% plot((vProduct))
% xlabel('Coefficient a_{i}')
% ylabel('log_{10} Product of a_{i}')
% hold off
% 
% figure_name = sprintf('%s : Sum of coefficient Scalars in C(f)',mfilename);
% figure('name',figure_name)
% hold on
% plot((vSum))
% xlabel('Coefficient a_{i}')
% ylabel('Sum of a_{i}')
% hold off

end
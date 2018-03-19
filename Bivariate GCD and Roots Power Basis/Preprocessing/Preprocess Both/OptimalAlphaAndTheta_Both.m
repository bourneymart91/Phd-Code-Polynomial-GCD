function [alpha,theta1,theta2] = OptimalAlphaAndTheta_Both(fxy_matrix, gxy_matrix,m,n)



% Define vector f for optimisation problem.
f = [1 -1 0 0 0];


%%
% get the degree of polynomial f and g
[m1, m2] = GetDegree_Bivariate(fxy_matrix);
[n1, n2] = GetDegree_Bivariate(gxy_matrix);

% Get the number of non-zeros in f(x,y) and g(x,y)
nNonZeros_fxy = GetNumNonZeros(m1,m2,m);
nNonZeros_gxy = GetNumNonZeros(n1,n2,n);

i1s_f = GetAsVector(diag(0:1:m1)*ones(m1+1,m2+1));
i2s_f = GetAsVector(ones(m1+1,m2+1)*diag(0:1:m2));
i1s_f = i1s_f(1:nNonZeros_fxy);
i2s_f = i2s_f(1:nNonZeros_fxy);


PartOne = ...
    [
        ones(nNonZeros_fxy,1) ...
        zeros(nNonZeros_fxy,1) ...
        -1 .* i1s_f(1:nNonZeros_fxy) ...
        -1 .* i2s_f(1:nNonZeros_fxy) ...
        zeros(nNonZeros_fxy,1)
    ];

i1s_g = GetAsVector(diag(0:1:n1) * ones(n1+1,n2+1));
i2s_g = GetAsVector(ones(n1+1,n2+1) * diag(0:1:n2));
i1s_g = i1s_g(1:nNonZeros_gxy);
i2s_g = i2s_g(1:nNonZeros_gxy);

PartTwo = ...
    [
        ones(nNonZeros_gxy,1) ...
        zeros(nNonZeros_gxy,1) ...
        -1 .*  i1s_g(1:nNonZeros_gxy) ...
        -1 .*  i2s_g(1:nNonZeros_gxy) ...
        -1 .* ones(nNonZeros_gxy,1)
    ];

PartThree = ...
    [
        zeros(nNonZeros_fxy,1) ...
        -1.*ones(nNonZeros_fxy,1) ...
        i1s_f(1:nNonZeros_fxy) ...
        i2s_f(1:nNonZeros_fxy) ...
        zeros(nNonZeros_fxy,1)
    ];

PartFour = ...
    [
        zeros(nNonZeros_gxy,1) ...
        -1.*ones(nNonZeros_gxy,1) ...
        i1s_g(1:nNonZeros_gxy) ...
        i2s_g(1:nNonZeros_gxy) ...
        ones(nNonZeros_gxy,1)
    ];



% Now build the vector b = [\lambda ; \mu ;  \rho ;  \tau]
fxy_vec     = GetAsVector(abs(fxy_matrix));
lambda_vec  = fxy_vec(1:nNonZeros_fxy); 
rho_vec     = fxy_vec(1:nNonZeros_fxy);

gxy_vec     = GetAsVector(abs(gxy_matrix));
mu_vec      = gxy_vec(1:nNonZeros_gxy);
tau_vec     = gxy_vec(1:nNonZeros_gxy);




% Get the index of the zero rows in lambda
index_zero_lambda = find(lambda_vec==0);
% Remove the corresponding zeros from lambda
lambda_vec(index_zero_lambda,:) = [];
% Remove the corresponding rows from PartOne Matrix
PartOne(index_zero_lambda,:) = [];


% Get the index of the zero rows in lambda
index_zero_mu = find(mu_vec==0);
% Remove the corresponding zeros from lambda
mu_vec(index_zero_mu,:) = [];
% Remove the corresponding rows from PartOne Matrix
PartTwo(index_zero_mu,:) = [];


% Get the index of the zero rows in rho
index_zero_rho = find(rho_vec==0);
% Remove the corresponding zeros from rho
rho_vec(index_zero_rho,:) = [];
% Remove the corresponding rows from PartOne Matrix
PartThree(index_zero_rho,:) = [];


% Get the index of the zero rows in tau
index_zero_tau = find(tau_vec==0);
% Remove the corresponding zeros from tau
tau_vec(index_zero_tau,:) = [];
% Remove the corresponding rows from PartOne Matrix
PartFour(index_zero_tau,:) = [];

A =-[PartOne; PartTwo; PartThree; PartFour];


b = -[log10(abs(lambda_vec)); log10(abs(mu_vec)); -log10(abs(rho_vec));-log10(abs(tau_vec))];


x = linprog(f,A,b);

try
    theta1 = 10^x(3);
    theta2 = 10^x(4);
    alpha = 10^x(5);
    %fprintf('Optimal theta 1 and theta 2 given by: \n  theta_{1}: %0.5e \n  theta_{2}: %0.5e',theta1,theta2)
catch
    fprintf('Failed to optimize\n')
    theta1 = 1;
    theta2 = 1;
    alpha = 1; 

end


end
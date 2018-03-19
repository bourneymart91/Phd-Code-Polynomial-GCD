function [theta1,theta2] = OptimalTheta(fxy_matrix, gxy_matrix)


% define vector f
f = [1 -1 0 0];

% get the degree of polynomial f and g
[m1, m2] = GetDegree_Bivariate(fxy_matrix);
[n1, n2] = GetDegree_Bivariate(gxy_matrix);

nEntries_f = (m1+1) * (m2+1);
nEntries_g = (n1+1) * (n2+1);

v_i1_f = GetAsVector(diag(0:1:m1) * ones(m1+1,m2+1));
v_i2_f = GetAsVector(ones(m1+1,m2+1) * diag(0:1:m2));

v_i1_g = GetAsVector(diag(0:1:n1) * ones(n1+1,n2+1));
v_i2_g = GetAsVector(ones(n1+1,n2+1) * diag(0:1:n2));


% Assemble the four submatrices of Matrix A

PartOne = ...
    [
    ones(nEntries_f,1) ,... 
    zeros(nEntries_f,1) ,...
    -1.* v_i1_f ,...
    -1.* v_i2_f ...
    ];

PartTwo = ...
    [
    ones(nEntries_g,1) ,...
    zeros(nEntries_g,1) ,...
    -1 .* v_i1_g ,...
    -1 .* v_i2_g ...
    ];

PartThree = ...
    [
    zeros(nEntries_f,1) ,...
    -1.* ones(nEntries_f,1) ,... 
    v_i1_f ,...
    v_i2_f ...
    ];

PartFour = ...
    [
    zeros(nEntries_g,1) ,...
    -1 .* ones(nEntries_g,1) ,...
    v_i1_g ,...
    v_i2_g ...
    ];

% Now build the vector b
lambda_vec = GetAsVector(fxy_matrix);
mu_vec = GetAsVector(gxy_matrix);
rho_vec = GetAsVector(fxy_matrix);
tau_vec = GetAsVector(gxy_matrix);


% % 
% % Removing any zero coefficients
% %
% %
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

A = [PartOne; PartTwo; PartThree; PartFour];


b = [log10(abs(lambda_vec)); log10(abs(mu_vec)); -log10(abs(rho_vec));-log10(abs(tau_vec))];


x = linprog(f,-A,-b);

try
    theta1 = 10^x(3);
    theta2 = 10^x(4);
    %fprintf('Optimal theta 1 and theta 2 given by: \n  theta_{1}: %0.5e \n  theta_{2}: %0.5e',theta1,theta2)
catch
    fprintf('Failed to optimize\n')
    theta1 = 1;
    theta2 = 1;


end


end
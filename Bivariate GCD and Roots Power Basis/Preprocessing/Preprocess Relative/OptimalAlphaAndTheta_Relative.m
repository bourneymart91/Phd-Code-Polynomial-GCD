function [alpha,theta1,theta2] = OptimalAlphaAndTheta_Relative(fxy, gxy)
% OptimalAlphaAndTheta(fxy_matrix, gxy_matrix)
%
% Obtain the optimal values of alpha, theta1 and theta2 for the Sylvester 
% matrix of the two polynomials f(x,y) and g(x,y)
% 
% Inputs.
% 
% fxy_matrix : Coefficients of the polynomial f(x,y)
%
% gxy_matrix : Coefficients of the polynomial g(x,y)


% Define vector f
f = [1 -1 0 0 0];

% get the degree of polynomial f and g

[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

nEntries_f = (m1 + 1) * (m2 + 1);
nEntries_g = (n1 + 1) * (n2 + 1);


v_i1_f = GetAsVector(diag(0:1:m1) * ones(m1+1,m2+1));
v_i2_f = GetAsVector(ones(m1+1,m2+1) * diag(0:1:m2));

v_i1_g = GetAsVector(diag(0:1:n1) * ones(n1+1,n2+1));
v_i2_g = GetAsVector(ones(n1+1,n2+1) * diag(0:1:n2));

PartOne = ...
    [
    ones(nEntries_f,1),...
    zeros(nEntries_f,1),...
    -1 .* v_i1_f ,...
    -1 .* v_i2_f ,...
    zeros(nEntries_f,1) ...
    ];

PartTwo = ...
    [
    ones(nEntries_g,1) ,...
    zeros(nEntries_g,1) ,...
    -1 .* v_i1_g ,...
    -1 .* v_i2_g ,...
    -1.* ones(nEntries_g,1) ...
    ];

PartThree = ...
    [
    zeros(nEntries_f,1) ,...
    -1 .* ones(nEntries_f,1) ,...
    v_i1_f ,...
    v_i2_f ,...
    zeros(nEntries_f,1) ...
    ];

PartFour = ...
    [
    zeros(nEntries_g,1) ,...
    -1 .* ones(nEntries_g,1) ,...
    v_i1_g ,...
    v_i2_g ,...
    ones(nEntries_g,1)
    ];


% Now build the vector b = [\lambda ; \mu ;  \rho ;  \tau]
lambda_vec = GetAsVector(abs(fxy));
mu_vec = GetAsVector(abs(gxy));
rho_vec = GetAsVector(abs(fxy));
tau_vec = GetAsVector(abs(gxy));


% 
% Get the index of the zero rows in lambda
index_zero_lambda = find(lambda_vec==0);

% Remove the corresponding zeros from lambda
lambda_vec(index_zero_lambda,:) = [];

% Remove the corresponding rows from PartOne Matrix
PartOne(index_zero_lambda,:) = [];


% Get the index of the zero rows in mu vector
index_zero_mu = find(mu_vec==0);

% Remove the corresponding zeros from mu vector
mu_vec(index_zero_mu,:) = [];

% Remove the corresponding rows from PartOne Matrix
PartTwo(index_zero_mu,:) = [];


% Get the index of the zero rows in rho
index_zero_rho = find(rho_vec==0);

% Remove the corresponding zeros from rho
rho_vec(index_zero_rho,:) = [];

% Remove the corresponding rows from PartOne Matrix
PartThree(index_zero_rho,:) = [];


% Get the index of the zero rows in tau vector
index_zero_tau = find(tau_vec==0);

% Remove the corresponding zeros from tau
tau_vec(index_zero_tau,:) = [];

% Remove the corresponding rows from PartOne Matrix
PartFour(index_zero_tau,:) = [];

A = -[PartOne; PartTwo; PartThree; PartFour];


b = -[log10((lambda_vec)); log10((mu_vec)); -log10((rho_vec));-log10((tau_vec))];

warning('off')
x = linprog(f,A,b);
warning('on')

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
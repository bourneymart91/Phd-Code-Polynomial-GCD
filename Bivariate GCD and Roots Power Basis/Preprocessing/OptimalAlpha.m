function [alpha] = OptimalAlpha(fxy_matrix, gxy_matrix)
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
%
% % Outputs
%
% alpha : Optimal value of \alpha, such that the two partitions of the
% Sylvester matrix S_{k}(f,\alpha*g) are balanced.

% Define vector f
f = [1 -1 0];

% Get the degree of polynomial f(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy_matrix);

% Get the degree of polynomial g(x,y) with respect to x and y.
[n1, n2] = GetDegree_Bivariate(gxy_matrix);

% Get the number of coefficients of f(x,y)
nEntries_f = (m1 + 1) * (m2 + 1);

% Get the number of coefficients of g(x,y)
nEntries_g = (n1 + 1) * (n2 + 1);

% Get the matrix whose columns are [i1 | i2] which are the powers of
% \theta_{1} and \theta_{2} corresponding to the vector of coefficients
% a_{i1,i2}.

% Get vector of i1
v_i1_f = GetAsVector(diag(0:1:m1) * ones(m1+1,m2+1));

% Get vector of i2
v_i2_f = GetAsVector(ones(m1+1,m2+1) * diag(0:1:m2));

% Get the vector i_{1} and i_{2} 
v_i1_g = GetAsVector(diag(0:1:n1) * ones(n1+1,n2+1));
v_i2_g = GetAsVector(ones(n1+1,n2+1) * diag(0:1:n2));

PartOne = ...
    [
    ones(nEntries_f,1),...
    zeros(nEntries_f,1),...
    zeros(nEntries_f,1) ...
    ];

PartTwo = ...
    [
    ones(nEntries_g,1) ,...
    zeros(nEntries_g,1) ,...
    -1.* ones(nEntries_g,1) ...
    ];

PartThree = ...
    [
    zeros(nEntries_f,1) ,...
    -1 .* ones(nEntries_f,1) ,...
    zeros(nEntries_f,1) ...
    ];

PartFour = ...
    [
    zeros(nEntries_g,1) ,...
    -1 .* ones(nEntries_g,1) ,...
    ones(nEntries_g,1)
    ];


% Now build the vector b = [\lambda ; \mu ;  \rho ;  \tau]
lambda_vec = GetAsVector(fxy_matrix);
mu_vec = GetAsVector(gxy_matrix);
rho_vec = GetAsVector(fxy_matrix);
tau_vec = GetAsVector(gxy_matrix);


% 
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


b =  [log10(abs(lambda_vec)); log10(abs(mu_vec)); -log10(abs(rho_vec));-log10(abs(tau_vec))];


x = linprog(f,-A,-b);

try
    alpha = 10^x(3);
catch
    fprintf([ mfilename 'Failed to optimize\n'])
    alpha = 1; 

end


end
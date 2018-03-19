function [alpha, theta1, theta2] = OptimalAlphaTheta(max_matrix_fxy, ...
    min_matrix_fxy, max_matrix_gxy, ...
    min_matrix_gxy, m, n)
%
%
% % Inputs
%
% max_matrix_fxy : (Matrix)
%
% min_matrix_fxy : (Matrix)
%
% max_matrix_gxy : (Matrix)
%
% min_matrix_gxy : (Matrix)
%
% m : (Int)
%
% n : (Int)
%
% % Outputs
%
% alpha : (Float)
%
% theta1 : (Float)
%
% theta2 : (Float)

% Obtain the optimal value of alpha and theta

% Define vector f
f = [1 -1 0 0 0];

% Assemble the four submatrices of Matrix A

% Get number of entries in f(x,y)
nEntries_fxy = nchoosek(m + 2, 2);

% Get number of entires in g(x,y)
nEntries_gxy = nchoosek(n + 2, 2);

v_i1_fxy = GetAsVector(diag(0 : 1 : m) * ones(m + 1, m + 1));
v_i2_fxy = GetAsVector(ones(m + 1, m + 1) * diag(0 : 1: m));
v_i1_fxy = v_i1_fxy(1 : nEntries_fxy);
v_i2_fxy = v_i2_fxy(1 : nEntries_fxy);

v_i1_gxy = GetAsVector(diag(0 : 1 : n) * ones(n + 1, n + 1));
v_i2_gxy = GetAsVector(ones(n + 1, n + 1) * diag(0 : 1 : n));
v_i1_gxy = v_i1_gxy(1 : nEntries_gxy);
v_i2_gxy = v_i2_gxy(1 : nEntries_gxy);


PartOne = ...
    [
    ones(nEntries_fxy,1) ...
    zeros(nEntries_fxy,1) ...
    -1 .* v_i1_fxy ...
    -1 .* v_i2_fxy ...
    zeros(nEntries_fxy,1)
    ];

PartTwo = ...
    [
    ones(nEntries_gxy,1) ...
    zeros(nEntries_gxy,1) ...
    -1.* v_i1_gxy ...
    -1.* v_i2_gxy ...
    -1.* ones(nEntries_gxy,1)
    ];

PartThree = ...
    [
    zeros(nEntries_fxy,1) ...
    -1.* ones(nEntries_fxy,1) ...
    v_i1_fxy ...
    v_i2_fxy ...
    zeros(nEntries_fxy,1)
    ];

PartFour = ...
    [
    zeros(nEntries_gxy,1) ...
    -1 .* ones(nEntries_gxy,1) ...
    v_i1_gxy ...
    v_i2_gxy ...
    ones(nEntries_gxy,1)
    ];



% Now build the vector b
lambda_vec = GetAsVector(max_matrix_fxy);
lambda_vec = lambda_vec(1:nEntries_fxy);

mu_vec = GetAsVector(max_matrix_gxy);
mu_vec = mu_vec(1:nEntries_gxy);

rho_vec = GetAsVector(min_matrix_fxy);
rho_vec = rho_vec(1:nEntries_fxy);

tau_vec = GetAsVector(min_matrix_gxy);
tau_vec = tau_vec(1:nEntries_gxy);

% % Find any zeros in the lambda vector
indeces = find(~lambda_vec);
PartOne(indeces,:) = [];
lambda_vec(indeces,:) = [];

% Find any zeros in the mu vector
indeces = find(~mu_vec);
PartTwo(indeces,:) = [];
mu_vec(indeces,:) = [];

% Find any zeros in the rho vector
indeces = find(~rho_vec);
PartThree(indeces,:) = [];
rho_vec(indeces,:) = [];

% Find any zeros in the tau vector
indeces = find(~tau_vec);
PartFour(indeces,:) = [];
tau_vec(indeces,:) = [];


b = [log10(lambda_vec); log10(mu_vec); -log10(rho_vec);-log10(tau_vec)]';

A = [PartOne; PartTwo; PartThree; PartFour];

warning('off')
x = linprog(f,-A,-b);
warning('on')

try
    alpha  = 10^x(5);
    theta1 = 10^x(3);
    theta2 = 10^x(4);
catch
    alpha = 1;
    theta1 = 1;
    theta2 = 1;
    
    
end

end
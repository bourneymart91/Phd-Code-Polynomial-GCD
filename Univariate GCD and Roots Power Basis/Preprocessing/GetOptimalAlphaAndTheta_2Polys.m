function [alpha, theta] = GetOptimalAlphaAndTheta_2Polys(fx,gx)

% Ensure that f(x) and g(x) are column vectors
if size(fx,2) >1  || size(gx,2)>2
    error('f(x) and g(x) must be column vectors')
end

f = [1 -1 0  0];

% Get degree of polynomial f
m = GetDegree(fx);

% Get degree of polynomial g
n = GetDegree(gx);

% Build the first partiion
Part1 = [ones(m+1,1) zeros(m+1,1) -1.*(0:1:m)' zeros(m+1,1)];

Part2 = [ones(n+1,1) zeros(n+1,1) -1.*(0:1:n)' -1.*ones(n+1,1)];

Part3 = [zeros(m+1,1) -1.*ones(m+1,1) (0:1:m)' zeros(m+1,1)];

Part4 = [zeros(n+1,1) -1.*ones(n+1,1) (0:1:n)' ones(n+1,1)];

lambda_vec = fx;
mu_vec = gx;
rho_vec = fx;
tau_vec = gx;



%% 
% Must first remove any zeros from lambda_vec, and corresponding rows in 
% Part one

% Get the index of the zero rows in lambda
index_zero_lambda = find(lambda_vec==0);
% Remove the corresponding zeros from lambda
lambda_vec(index_zero_lambda,:) = [];
% Remove the corresponding rows from PartOne Matrix
Part1(index_zero_lambda,:) = [];

% Get the index of the zero rows in lambda
index_zero_mu = find(mu_vec==0);
% Remove the corresponding zeros from lambda
mu_vec(index_zero_mu,:) = [];
% Remove the corresponding rows from PartOne Matrix
Part2(index_zero_mu,:) = [];

% Get the index of the zero rows in lambda
index_zero_rho = find(rho_vec==0);
% Remove the corresponding zeros from lambda
rho_vec(index_zero_rho,:) = [];
% Remove the corresponding rows from PartOne Matrix
Part3(index_zero_rho,:) = [];

% Get the index of the zero rows in lambda
index_zero_tau = find(tau_vec==0);
% Remove the corresponding zeros from lambda
tau_vec(index_zero_tau,:) = [];
% Remove the corresponding rows from PartOne Matrix
Part4(index_zero_tau,:) = [];



A = -[Part1; Part2; Part3; Part4];


b = -[log10(abs(lambda_vec)); log10(abs(mu_vec)); -log10(abs(rho_vec));-log10(abs(tau_vec))];

warning('off')
x = linprog(f,A,b);
warning('on')


try
    theta = 10^x(3);
    alpha = 10^x(4);
catch
    fprintf('Failed to optimize\n')
    theta = 1;
    alpha = 1; 

end

end
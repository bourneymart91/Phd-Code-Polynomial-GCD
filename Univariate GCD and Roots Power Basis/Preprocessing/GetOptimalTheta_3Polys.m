function [theta] = GetOptimalTheta_3Polys(fx, gx, hx)
%
% % Inputs
%
% [fx, gx, hx] : Coefficients of polynomials f(x), g(x) and h(x)

% Ensure that f(x) and g(x) are column vectors
if size(fx,2) >1  || size(gx,2)>2
    error('f(x) and g(x) must be column vectors')
end

f = [1 -1 0];

% Get degree of polynomial f(x)
m = GetDegree(fx);

% Get degree of polynomial g(x)
n = GetDegree(gx);

% Get degree of polynomial h(x)
o = GetDegree(hx);

% Build the first partiion
Part1 = [ones(m+1,1) zeros(m+1,1) -1.*(0:1:m)' ];

Part2 = [ones(n+1,1) zeros(n+1,1) -1.*(0:1:n)' ];

Part3 = [ones(o+1,1) zeros(o+1,1) -1.*(0:1:o)' ];

Part4 = [zeros(m+1,1) -1.*ones(m+1,1) (0:1:m)' ];

Part5 = [zeros(n+1,1) -1.*ones(n+1,1) (0:1:n)' ];

Part6 = [zeros(o+1,1) -1.*ones(o+1,1) (0:1:o)' ];


f_vec_a = fx;
g_vec_a = gx;
h_vec_a = hx;

f_vec_b = fx;
g_vec_b = gx;
h_vec_b = hx;


%% 
% Must first remove any zeros from lambda_vec, and corresponding rows in 
% Part one

% Get the index of the zero rows in lambda
index_zero_f_vec_a = find(f_vec_a==0);
% Remove the corresponding zeros from lambda
f_vec_a(index_zero_f_vec_a,:) = [];
% Remove the corresponding rows from PartOne Matrix
Part1(index_zero_f_vec_a,:) = [];

% Get the index of the zero rows in lambda
index_zero_g_vec_a = find(g_vec_a==0);
% Remove the corresponding zeros from lambda
g_vec_a(index_zero_g_vec_a,:) = [];
% Remove the corresponding rows from PartOne Matrix
Part2(index_zero_g_vec_a,:) = [];

% Get the index of the zero rows in lambda
index_zero_h_vec_a = find(h_vec_a==0);
% Remove the corresponding zeros from lambda
g_vec_a(index_zero_h_vec_a,:) = [];
% Remove the corresponding rows from PartOne Matrix
Part3(index_zero_h_vec_a,:) = [];




% Get the index of the zero rows in lambda
index_zero_f_vec_b = find(f_vec_b==0);
% Remove the corresponding zeros from lambda
f_vec_b(index_zero_f_vec_b,:) = [];
% Remove the corresponding rows from PartOne Matrix
Part4(index_zero_f_vec_b,:) = [];

% Get the index of the zero rows in lambda
index_zero_g_vec_b = find(g_vec_b==0);
% Remove the corresponding zeros from lambda
g_vec_b(index_zero_g_vec_b,:) = [];
% Remove the corresponding rows from PartOne Matrix
Part5(index_zero_g_vec_b,:) = [];


% Get the index of the zero rows in lambda
index_zero_h_vec_b = find(h_vec_b==0);
% Remove the corresponding zeros from lambda
g_vec_b(index_zero_h_vec_b,:) = [];
% Remove the corresponding rows from PartOne Matrix
Part6(index_zero_h_vec_b,:) = [];


A = -[Part1; Part2; Part3; Part4; Part5; Part6];


b = -[log10(abs(f_vec_a)); log10(abs(g_vec_a)); log10(abs(h_vec_a)); -log10(abs(f_vec_b)); -log10(abs(g_vec_b)); -log10(abs(h_vec_b))];

warning('off')
x = linprog(f,A,b);
warning('on')


try
    theta = 10^x(3);
catch
    fprintf('Failed to optimize\n')
    theta = 1;
    alpha = 1; 

end

end
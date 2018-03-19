function theta = GetOptimalTheta(arr_fx)
% Get optimal value of theta for the matrix.
%
%
% Inputs
%
% arr_fx : (Array of Vectors) Each cell in the array contains coefficients
% of the polynomial f_{i}(x)
%
% Outputs
%
% theta : (Float) Optimal value of \theta



% Get number of polynomials in the array
nPolynomials_arr_fx = size(arr_fx, 1);


% Get vector of degrees of polynomials f_{i}(x)
vDegree_arr_fx = zeros(nPolynomials_arr_fx,1);

for i = 1 : 1 : nPolynomials_arr_fx
    vDegree_arr_fx(i) = GetDegree(arr_fx{i});
end

% Get vector of degrees of polynomials h_{i}(x)
vDeg_arr_hx = vDegree_arr_fx(1:end-1) - vDegree_arr_fx(2:end);

% Initialise arrays to store maximum and minimum entry of each coefficient
% of each polynomial
arrF_max = cell(nPolynomials_arr_fx, 1);
arrF_min = cell(nPolynomials_arr_fx, 1);

% For each polynomial f_{1},...,f_{d}, note we exclude f_{0} from this,
% since f_{0} does not appear in the LHS matrix.
for i = 2 : 1 : nPolynomials_arr_fx
    
    % Get polynomial f_{i} from the set g containing all f_{i}
    fx = arr_fx{i};
    
    % Assign empty vectors for the max and minimum values of each
    % coefficient in F.
    
    
    n = vDeg_arr_hx(i-1);
    
    [vF_max, vF_min] = GetMaxMin(fx, n);
    
    arrF_max{i} = vF_max;
    arrF_min{i} = vF_min;
    
end

theta = MinMaxOpt(arrF_max, arrF_min);

end

function theta = MinMaxOpt(arrF_max, arrF_min)
%
% This function computes the optimal value theta for the preprocessing
% opertation as part of block deconvolution
%
% F_max   :  (Array of Vectors)
%
% F_min   :  A vector of length m+1, such that F_min(i) stores the
%            element of minimum magnitude of S(f,g)Q that contains the
%            coefficient a(i) of f, i=1,...,m+1.
%

% Get the number of polynomials
nPolynomials_arr_fxy = size(arrF_max,1);



f = [1 -1 0];

Ai = cell(nPolynomials_arr_fxy, 1);
Bi = cell(nPolynomials_arr_fxy, 1);


for i = 1 : 1 : nPolynomials_arr_fxy
    
    % Get the max of each coefficient of polynomial fw
    fw_max = arrF_max{i,1};
    
    % Get Degree of the polynomial
    m = GetDegree(fw_max);
    
    % Build the matrix A_{i}
    Ai{i,1} = [ones(m + 1, 1) zeros(m + 1, 1)   -(0 : 1 : m)'];
    
    % Build the matrix B_{i}
    Bi{i,1} = [zeros(m + 1, 1) -ones(m + 1, 1) (0 : 1 : m)'];
    
end

Part1 = cell2mat(Ai);

Part2 = cell2mat(Bi);

A = ...
    [
    Part1;
    Part2
    ];

% ------------------------------------------------------------------------

% Get the array of entries F_max_{i} as a vector
v_F_max = cell2mat(arrF_max);
v_F_min = cell2mat(arrF_min);

b = [log10(v_F_max); -log10(v_F_min)];



%--------------------------------------------------------------------------

% Solve the linear programming problem and extract alpha and theta
% from the solution vector x.
try
    
    x = linprog(f, -A, -b);
    theta = 10^x(3);
    fprintf([mfilename ' : ' sprintf('Optimal theta : %2.4f \n', theta)])
    
catch
    
    fprintf('Error Calculating Optimal value of theta\n');
    theta = 1;
    
end

end


function [vF_max, vF_min] = GetMaxMin(fx,n)
%
% % Inputs
%
% fx : (Vector) Coefficients of f(x)
%
% n : (Int) Degree of h(x)
%
% % Outputs
%
% vF_Max : (Vector) Contains maximum occurence of each coefficient a_{i} in
% the matrix T_{n}(f)
%
% vF_Min : (Vector) Contains minimum occurence of each coefficient a_{i}
% in the matrix T_{n}(f)

% Get degree
m = GetDegree(fx);

vF_max = zeros(m + 1, 1);
vF_min = zeros(m + 1, 1);

% For each coefficient a_{j} in f_{i+1}
for i = 0 : 1 : m
    
    % Get the coefficient a_{j} of polynomial f_{i}
    aij = fx(i+1);
    
    % initialise a vector to store all the a_{i} for j = 0,...,n
    vec_x = zeros(n + 1, 1);
    
    % For each occurence of the coefficient ai_j in the columns of C_{n_{i}}(f_{i})
    for j = 0 : 1 : n
        
        vec_x(j+1) = aij .* nchoosek(m, i) * nchoosek(n, j) ./ nchoosek(m + n, i + j);
        
    end
    
    % Get max entry of each coefficient.
    vF_max(i+1) = max(abs(vec_x));
    
    % Get min entry of each coefficient.
    vF_min(i+1) = min(abs(vec_x));
    
end

end
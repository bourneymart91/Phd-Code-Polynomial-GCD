function [th1, th2] = GetOptimalTheta(arr_fxy)
%
% % Inputs
%
% arr_fxy : (Array of Matrices) Array of matrices containing coefficients
% of each of the polynomials f_{i}(x,y)
%
% % Outputs
%
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}



% GetOptimalTheta(arr_fxy)
f = [1 -1 0 0];

% Get the number of polynomials in the array
nPolynomials_arr_fxy = size(arr_fxy,1);

% Get the degree of the array of polynomials
vDeg_arr_fxy = zeros(nPolynomials_arr_fxy,1);

for i = 1:1:nPolynomials_arr_fxy
    vDeg_arr_fxy(i) = GetDegree_Bivariate(arr_fxy{i});
end

% Get degree of array of polynomials h_{i}(x,y)
vDeg_arr_hxy = vDeg_arr_fxy(1:end-1) - vDeg_arr_fxy(2:end);




%--------------------------------------------------------------------------

% Create an array which stores the maximum of each coefficient from each
% polynomial f_{i}(x,y)
arr_max_fxy = cell(nPolynomials_arr_fxy,1);
arr_min_fxy = cell(nPolynomials_arr_fxy,1);

for i = 2 : 1 : nPolynomials_arr_fxy
    
    % For each polynomial in the array , get the maximum and minimum entry
    % of each coefficient
    
    [arr_max_fxy{i,1} , arr_min_fxy{i,1} ] = ...
        GetMaxMin_fxy(arr_fxy{i},vDeg_arr_hxy(i-1));
    
end

% Initialise cell arrays rho and tau
arrRho = cell(nPolynomials_arr_fxy,1);
arrTau = cell(nPolynomials_arr_fxy,1);

for i = 2 : 1 : nPolynomials_arr_fxy
    
    m = GetDegree_Bivariate(arr_fxy{i});
    nCoefficients_fxy = nchoosek(m+2,2);
    
    temp_max = GetAsVector(arr_max_fxy{i});
    temp_max = temp_max(1:nCoefficients_fxy);
    
    temp_min = GetAsVector(arr_min_fxy{i});
    temp_min = temp_min(1:nCoefficients_fxy);
    
    arrRho{i} = temp_max;
    arrTau{i} = temp_min;
    
end

rho_vec = cell2mat(arrRho);
tau_vec = cell2mat(arrTau);


% Build right hand side vector
RHS_Vec = [rho_vec ; - tau_vec];





% -------------------------------------------------------------------------

% Code for constructing LHS matrix


% For each of the polynomials f_{1}...,f_{d}
arr_A = cell(nPolynomials_arr_fxy,1);
arr_B = cell(nPolynomials_arr_fxy,1);

for i = 2:1:nPolynomials_arr_fxy
    
    % Get degree of f_{i}(x)
    m = vDeg_arr_fxy(i);
    
    % Get number of coefficients in f_{i}(x)
    nCoefficients_fxy = nchoosek(m + 2, 2);
    
    % Build part of the matrix
    arr_A{i,1} = [ones(nCoefficients_fxy,1) zeros(nCoefficients_fxy,1) -1.*GetPairs(m) ];
    arr_B{i,1} = [zeros(nCoefficients_fxy,1) -1.*ones(nCoefficients_fxy,1) 1.*GetPairs(m)];
    
end


% Build left hand side matrix
A = cell2mat(arr_A);
B = cell2mat(arr_B);
LHS_Matrix = [A; B];

% ------------------------------------------------------------------------



x = linprog(f, -LHS_Matrix, -RHS_Vec);


try
    th1 = 10^x(3);
    th2 = 10^x(4);
    
    fprintf('Optimal theta_{1} : %f \n', th1)
    fprintf('Optimal theta_{2} : %f \n', th1)
    
    
    
catch
    fprintf('Failed to optimize\n')
    th1 = 1;
    th2 = 1;
end

end












function vPairs = GetPairs(m)
%
% % Inputs
%
% m : (Int) Degree of polynomial f(x,y)
%
% % Outputs
%
% vPairs : (Matrix) 2 column matrix containing [i_{1}, i_{2}] pairs.

% Get number of coefficients in f(x,y)
nCoefficients_fxy = nchoosek(m + 2, 2);

% Get index for i1
mat_i1 = diag(0:1:m) * ones(m + 1, m + 1);
v_i1 = GetAsVector(mat_i1);
v_i1 = v_i1(1 : nCoefficients_fxy);

% Get index for i2
mat_i2 = ones(m + 1, m + 1) * diag(0 : 1 : m);
v_i2 = GetAsVector(mat_i2);
v_i2 = v_i2(1:nCoefficients_fxy);

% Create set of pairs
vPairs = [v_i1 v_i2];

end

function [max_matrix,min_matrix] = GetMaxMin_fxy(fxy, n)
% Get maximum occurence of each coefficient
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% n : (Int) Total degree of polynomial g(x,y)
%
% % Outputs
%
% max_matrix : (Matrix) Max occurence of each coefficient in T_{n}(f)
%
% min_matrix : (Matrix) Min occurence of each coefficient in T_{n}(f)


% Get the degree of f(x,y)
m = GetDegree_Bivariate(fxy);

% Initialise matrices
max_matrix = zeros(m + 1, m + 1);
min_matrix = zeros(m + 1, m + 1);

for i = 0:1:m
    
    for i1 = i:-1:0
        
        i2 = i - i1;
        
        % Get the coefficient
        a_i1i2 = fxy(i1 + 1, i2 + 1);
        
        numerator_binomial_1 = Trinomial(m, i1, i2);
        
        % Initialise matrix to store each occurence of a_{i1,i2}
        ai1i2 = zeros(n + 1, n + 1);
        
        for j = 0:1:n
            for j1 = j:-1:0
                
                j2 = j - j1;
                
                numerator_binomial_2 = Trinomial(n, j1, j2);
                
                denominator_binomial_1 = Trinomial(m + n, i1 + j1, i2 + j2);
                
                ai1i2(j1 + 1,j2 + 1) = ...
                    abs(a_i1i2 * numerator_binomial_1 * numerator_binomial_2 ./ denominator_binomial_1);
                
            end
        end
        
        
        % Get Vector of each occurence of a_{i1,i2}
        nCoefficients_hxy = nchoosek(n + 2, 2);
        v_ai1i2 = GetAsVector(ai1i2);
        v_ai1i2 = v_ai1i2(1:nCoefficients_hxy);
        
        max_ai1i2 = max(v_ai1i2);
        min_ai1i2 = min(v_ai1i2);
        
        max_matrix(i1 + 1, i2 + 1) = max_ai1i2;
        min_matrix(i1 + 1, i2 + 1) = min_ai1i2;
        
    end
end






end
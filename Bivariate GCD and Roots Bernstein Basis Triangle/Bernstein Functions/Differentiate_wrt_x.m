function [dfxy] = Differentiate_wrt_x(fxy)
% Differentiate_wrt_x(fxy)
%
% Differentiate with respect to x, the bivariate polynomial f(x,y) which is
% given in Bernstein form.
% 
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y).
%
% % Outputs
%
% dfxy : (Matrix) The coefficients of the polynomial \frac{\partial f}{\partial x}

% Get the degree of f(x,y)
[m1,~] = GetDegree_Bivariate(fxy);
m = m1;

% Get number of coefficients of f(x,y)
nCoefficients_fxy = nchoosek(m+2,2);

% Get vector of coefficients of polynomial f(x,y)
v_fxy = GetAsVector(fxy);

% Remove zeros from end of vector to get coefficients of f(x,y)
v_fxy = v_fxy(1 : nCoefficients_fxy);

% Initialise a temp matrix, consisting of zeros
temp_mat = zeros(m+1, m+1);

% Get the number of coefficients in the derivative
nCoefficients_dfxy = nchoosek(m+1, 2);

% Initialise a matrix which transforms the coefficients of the polynomial
% f(x,y) to the coefficients of the derivative of f(x,y).
trans_mat = zeros(nCoefficients_dfxy, nCoefficients_fxy);

% Initialise a count
count = 1;


for k = 0:1:m-1
    for i = k:-1:0
    
        mat = temp_mat;
        
        j = k-i;
    
        mat(i+1,j+1) = -m;
        mat(i+2,j+1) = m;
    
        temp_vec = GetAsVector(mat)';
        temp_vec = temp_vec(1:nCoefficients_fxy);
        
        trans_mat(count,:) = temp_vec;
        count = count + 1;
    end
end

% Get nonzero Coefficients of f(x,y)
v_dfxy = trans_mat * v_fxy;

if m > 1
    nZeros = nchoosek(m,2);
else
    nZeros = 0;
end

% Append zeros to make dfxy fit into a square matrix of size (m-1) x (m-1)
v_dfxy = ...
    [
    v_dfxy;
    zeros(nZeros,1);
    ];


dfxy = GetAsMatrix(v_dfxy,m-1,m-1);


end
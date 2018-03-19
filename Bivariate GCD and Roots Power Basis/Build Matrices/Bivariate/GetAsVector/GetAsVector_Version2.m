function [v_fxy] = GetAsVector_NewMethod(fxy)
% Given the matrix of coefficients of polynomial f(x,y), insert the
% coefficients into a vector reading down the columns from left to right.
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% % Outputs
%
% v_fxy : (Vector) Vector of coefficients of f(x,y)

% Get degree of f(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);

% Get number of coefficients in f(x,y)
nCoefficients = (m1+1) * (m2+1);

% Initialise the vector to store coefficients of f(x,y)
v_fxy = zeros(nCoefficients, 1);

% For each column
for j = 1:1:(m2+1)
   
    % Get column of coefficients
    temp_vec = fxy(:,j);
    
    % Get start and end point of vector to insert the column of
    % coefficients into.
    start_index = ((j-1) * (m1+1)) + 1;
    end_index = (j) * (m1+1);
    
    % Insert coefficients from jth column into the vector
    v_fxy(start_index : end_index) = temp_vec; 
    
    
end



end
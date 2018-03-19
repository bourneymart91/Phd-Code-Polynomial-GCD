function [max_matrix, min_matrix] = GetMaxMin(fxy, m, n_k)
% Get the maximum and minimum of each entry of f(x,y) in the Sylvester
% matrix.
%
% Inputs
%
% fxy : The Coefficients of polynomial f(x,y) in standard bernstein basis.
%       Given in matrix form so that the rows are in terms of x basis
%       elements and the columns are y basis elements.
%
% m : Total degree of f(x,y)
%
% n_k : Total degree of v_{k}(x,y)
%
%
% % Outputs
%
% max_matrix : (Matrix) Matrix containing maximum occurence of each
% coefficient of f(x,y) in the subresultant 
%
% min_matrix : (Matrix) Matrix containing minimum occurence of each
% coefficietn of f(x,y) in the subresultant.
%
%
%

% Take the absolute values of the coefficients
fxy = abs(fxy);

% Build a matrix which stores the maximum values of each coefficient
max_matrix = zeros(m + 1, m + 1);
 
% Build a matrix which stores the minimum values of each coefficient
min_matrix = zeros(m + 1, m + 1);


% for each coefficient a_{i1,i2} - Get the maximum and minimum entry of
% a_{i1,i2} in DTQ


% For each diagonal of the matrix containing coefficients a_{i1,i2}
for k = 0 : 1: m
    
    % Get the index i_{1} for the row of the matrix
    for i1 = k : -1 : 0
        
        % Get the index i_{2} for the col of the matrix
        i2 = k - i1;
                
        % Get the coefficient
        ai1i2 = fxy(i1 + 1, i2 + 1);
        
        % Get max and min entry of the coefficient
        [max_ai1i2, min_ai1i2] = GetMaxMin2(ai1i2, i1, i2, m, n_k);
        
        max_matrix(i1+1,i2+1) = max_ai1i2;
        min_matrix(i1+1,i2+1) = min_ai1i2;
        
    end
end


end

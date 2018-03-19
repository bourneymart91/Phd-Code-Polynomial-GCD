
function f_vec = GetAsVector_Version1(fxy_matrix)
% Given the matrix of coefficients of the polynomial f(x,y), return the 
% vector of the coefficients such that the order is increasing and for 
% coefficients of the same order the highest power of x is first.
%
% % Inputs
%
% fxy_matrix : (Matrix) Contains coefficients of polynomial f(x,y)
%
% % Outputs
%
% f_vec : (Vector) Contains coefficients of polynomial f(x,y)

% Get degree of polynomial f(x,y) with respect to x
[m1, m2] = GetDegree_Bivariate(fxy_matrix);

% Initialise a count
count = 1;

% Initialise an empty vector for coefficients of f(x,y)
f_vec = zeros((m1+1)*(m2+1),1);

% Get the number of diagonals in the matrix of coefficients of f(x,y)
nDiagonals = (m1+1)+(m2+1)-1;

% For each diagonal of the matrix
for tot = 0:1:nDiagonals
    for i = tot:-1:0
        j = tot - i;
        
        if(i > m1 || j > m2)
            
        else
            f_vec(count) = fxy_matrix(i+1,j+1);
            count = count + 1;
        end
        
    end
    
end


end
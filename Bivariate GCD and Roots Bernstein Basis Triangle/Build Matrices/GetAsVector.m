function f_vec = GetAsVector(fxy_matrix)
% GetAsVector(fxy_matrix)
%
% Given the polynomial f in the bernstein basis, whose coefficients are in
% matrix form, obtain the vector of the coefficients such that the order is
% increasing and the higher power of x is first.
%
% % Inputs.
%
% fxy_matrix : (Matrix) Coefficients of polynomial f(x,y)
%
% % Outputs.
% 
% f_vec : (Vector) Coefficients of polynomial f(x,y)

% Get degree of polynomial f(x,y) with respect to x
[m1, m2] = GetDegree_Bivariate(fxy_matrix);

% Initialise a count
count = 1;

% Initialise an empty vector for coefficients of f(x,y)
f_vec = zeros((m1 + 1) * (m2 + 1), 1);

% Get the number of diagonals in the matrix of coefficients of f(x,y)
nDiagonals = (m1 + 1) + (m2 + 1) -1;


for tot = 0 : 1 : nDiagonals
    for i = tot : -1 : 0
        j = tot - i;
        
        if( i > m1 || j > m2)
        else
            f_vec(count) = fxy_matrix(i + 1, j + 1);
            count = count + 1;
        end
        
    end
    
end


end
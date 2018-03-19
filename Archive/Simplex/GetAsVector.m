function f_vec = GetAsVector(fxy_matrix)
% Given the polynomial f in the bernstein basis, whose coefficients are in
% matrix form, obtain the vector of the coefficients such that the order is
% increasing and the higher power of x is first.

%% Get Structure of f(x,y)

% Get degree of polynomial f(x,y) with respect to x
[m1,m2] = GetDegree(fxy_matrix);

% Initialise a count
count = 1;

% Initialise an empty vector for coefficients of f(x,y)
f_vec = zeros((m1+1)*(m2+1),1);

% Get the number of diagonals in the matrix of coefficients of f(x,y)
num_diags = (m1+1)+(m2+1)-1;


for tot = 0:1:num_diags
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
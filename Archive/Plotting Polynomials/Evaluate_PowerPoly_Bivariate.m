
function [sum] = Evaluate_PowerPoly_Bivariate(t1,t2,fxy)
% % Given a bivariate polynomial in matrix form, evaluate the polynomial
% for a given point (t1,t2)

% % Inputs

% t1 :  Evaluation point in x

% t2 :  Evaluation point in y

% fxy:  Matrix of polynomial coefficients.

% % Outputs

% sum:  Value of the surface function at point (t1,t2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the degree of the polynomial with respect to x
m1 = size(fxy,1)-1;

% Get the degree of the polynomial with respect to y
m2 = size(fxy,2)-1;

% Perform a summing function to evaluate the curve at the given point
sum = 0;
% for each row, i = 0,..,m_{1}
for i = 0:1:m1
    % for each column j = 0,..,m_{2}
    for j = 0:1:m2
        temp_val = fxy(i+1,j+1)* (t1^(m1-i))  * (t2^(m2-j));
        sum = sum + temp_val;
    end
    
end

end

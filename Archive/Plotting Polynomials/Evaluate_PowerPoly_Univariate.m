function [sum] = Evaluate_PowerPoly_Univariate(t1,fx)
% Given the function fx, evaluate at point t1

% % Inputs 

% t1:   Evaluation point

% fx:   Vector of coefficients of fx where leading coefficient is of
%       highest degree. a_{i}y^{m-i}

% % Outputs

% sum

% Get the degree of the input polynomial fx
m = length(fx) -1;

% Initialise a sum
sum = 0;

% for each entry in the coefficient vector
for i = 0:1:m
    sum = sum + (fx(i+1) * (t1^(m-i)));
end

end
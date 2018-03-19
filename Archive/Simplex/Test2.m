function [] = Test2

% Initialise symbols
syms b_0 b_1 b_2 b_3 b_4

%
% Get example
degree = 4;

switch degree
    case 1
        fx = [b_0; b_1];
        
    case 2
        fx = [b_0; b_1; b_2];
    case 3
        fx = [b_0; b_1; b_2; b_3]
    case 4
        fx = [b_0; b_1; b_2; b_3; b_4]
end

% Get degree
m = GetDegree(fx);

% Initialise symbolic x
x = sym('x');

% %
% %
% Get symbolic expression 
sum = 0;
for k = 0:1:m
    
    sum = sum + ...
        (...
        fx(k+1) * nchoosek(m,k) *...
        x^(k) * (1-x)^(m-k)...
        );
    
end

% % 
% %

my_exp = expand(sum)
collect(my_exp,[x])
end
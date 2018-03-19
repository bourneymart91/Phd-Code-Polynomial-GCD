function [fx] = GetCoefficients(root_mutl_arr_f)
% Given the roots and multiplicities of a polynomial f(x), calculate the
% coefficients
%
%   Inputs.
%
%   roots_f : Matrix consisting of rows of [root multiplicity] pairs
%

% Get the number of distinct roots
[nDistinctRoots,~] = size(root_mutl_arr_f);

% Initialise the polynomials
Poly = 1;

% for each distinct root
for i = 1:1:nDistinctRoots
    
    % Get the root
    r = root_mutl_arr_f(i,1);
    
    % Get the multiplicity
    m = root_mutl_arr_f(i,2);
    
    % Multiply poly by the factor (x-r) m times
    for j = 1:1:m
        temp = [1;-r];
        Poly = conv(Poly,temp);
    end
    
end

fx = flipud(Poly);
end
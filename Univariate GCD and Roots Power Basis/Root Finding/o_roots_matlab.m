function [root_multiplicity_matrix] = o_roots_matlab(fx)
% Given a vector of coefficients of a polynomail in Bernstein form. Return
% the set of roots given by zheng MultRoots() function.
%
% % Inputs.
%
% fx : Column vector of coefficients of the polynomial f(x) 
%
% % Outputs.
%
% roots_mult_array : Roots and multiplicities of f(x) as calculated by 
% the zheng MultRoots function.
%
%
%


% Build the vector of corresponding binomial coefficients
roots_calc = roots(flipud(fx));

% Get the number of roots
[nRoots,~] = size(roots_calc);

root_multiplicity_matrix = [roots_calc(:,1) ones(nRoots,1)];





end


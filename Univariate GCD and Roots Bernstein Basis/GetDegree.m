function m = GetDegree(fx)
% GetDegree(m)
%
% Get the degree of the polynomial f(x), whose coefficients are given as a
% vector.
%
% Inputs.
%
% fx : (Vector) Coefficients of vector f(x).
%
% % Outputs
%
% m : (Int) Degree of polynomial f(x)



% Get the degree.
m = size(fx, 1) - 1;

end
function [roots_Bb] = o_roots_multroot(fx)
% Given a vector of coefficients of a polynomail in Bernstein form. Return
% the set of roots given by zheng MultRoots() function.
%
% % Inputs.
%
% fx : (Vector) Vector of coefficients of the polynomial f(x) in Bernstein
%   form.
%
% % Outputs.
%
% roots_Bb : Roots of f(x) as calculated by the zheng MultRoots function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %

addpath 'Root Finding Methods/multroot/multroot'

% Get the degree of polynomial f(x) 
m = GetDegree(fx);

% Build the vector of corresponding binomial coefficients
Bi_m = GetBinomials(m);

fx_bi = fx.* Bi_m;

fx_bi = flipud(fx_bi);

roots_calc = multroot(fx_bi);

% convert roots wrt t, Bernstein basis
roots_Bb = [roots_calc(:,1)./(1+roots_calc(:,1)) roots_calc(:,2)];

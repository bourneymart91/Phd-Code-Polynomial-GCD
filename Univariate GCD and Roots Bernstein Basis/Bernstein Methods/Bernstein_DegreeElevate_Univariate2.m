function [deg_elv_f] = DegreeElevate_Univariate2(fx,r)
% Function performs degree elevation of a polynomial f(x) in the Bernstein
% basis.
%
% % Input.
%
% fx : The coefficients of polynomial f(x) of degree m. 
%
% r : Number of degree elevations such that the output polynomial is of 
% degree m + r.

% Get degree of polynomial f
m = GetDegree(fx);

% To degree elevate we multiply f(x) by a polynomial g(x), where the
% coefficients of g(x) are all 1.

gx = ones(r+1,1);

T1 = BuildT1(gx,m);
D = BuildD(m,r);
Q = BuildQ1(m);

DTQ = D*T1*Q;

deg_elv_f = DTQ * fx ;




end
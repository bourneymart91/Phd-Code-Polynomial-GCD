function [f_root_mult_array,g_root_mult_array] = BuildRandomPolynomials(m,n,t,intvl_low,intvl_high)
% Given the degree of two polynomials, construct two sets of roots
% 
% Inputs.
% 
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% t : Degree of polynomial d(x)
%
% intvl_low : Lowest allowed root value
%
% intvl_high : Highest allowed root value
%
% Outputs.
%
%
% f_root_mult_arr : root and multiplicity array for polynomial f.
%
% g_root_mult_arr : root and multiplicity array for polynomial g.

a = intvl_low;
b = intvl_high;

format long

% % 
% % 
% %

% Get the multiplicity structure of the roots of the GCD d of degree
% t. t = t1 + t2 + ... + t_{r}
% We want more lower multiplicity roots. so skew this way.
% 
mult_arr_d = GetMultiplicities(t);

% Get the number of roots of d
nDistinctRoots_d = length(mult_arr_d);

% %
% %
% %
% Get multiplicity structure for the roots of f(x,y) excluding those of d(x,y).
% The remaining roots of f are the roots of u(x,y)
mult_arr_u = GetMultiplicities(m-t);


if m-t ~=0
    nDistinctRoots_u = length(mult_arr_u);
    mult_arr_f = [mult_arr_d mult_arr_u];
else
    nDistinctRoots_u = 0;
    mult_arr_f = mult_arr_d;
end


% % 
% %
% %
% % Get multiplicity structure of g
% initialise a probability vector so that lower multiplicities are
% preferred

mult_arr_v = GetMultiplicities(n-t);

if (n-t ~= 0)
    nDistinctRoots_v = length(mult_arr_v);
    mult_arr_g = [mult_arr_d mult_arr_v];
else
    nDistinctRoots_v = 0;
    mult_arr_g = mult_arr_d;
end

% Get a set of unique roots
% the 1000 and 1000 contain the roots to the unit interval
detail = 100;
format 'long';



roots = a + randperm(detail,nDistinctRoots_d + nDistinctRoots_v + nDistinctRoots_u)./(detail./(b-a));


roots_d = roots(1:nDistinctRoots_d);
roots(1:nDistinctRoots_d) = [];
roots_f = roots(1:nDistinctRoots_u);
roots(1:nDistinctRoots_u) = [];
roots_g = roots(1:nDistinctRoots_v);


f_root_mult_array = [[roots_d'; roots_f'] mult_arr_f'];
g_root_mult_array = [[roots_d'; roots_g'] mult_arr_g'];


end


function [mult_arr] = GetMultiplicities(t)

global SETTINGS

mult_arr = [];

% Get a probability distribution.
prob_arr = zeros(1,t);
for i = 1:1:t
    prob_arr(i) = i./ nchoosek(t+1,2);
end
prob_arr = fliplr(prob_arr);
rng(SETTINGS.SEED);


% Get the multiplicity structure of d(x,y)
total = 0;
i = 1;
while total < t
    r = rand;
    prob = prob_arr;
    x = sum(r >= cumsum([0, prob]));
    if (total + x) <= t
        mult_arr(i) = x;
        total = total + x;
        i = i+1;
    end
end
end
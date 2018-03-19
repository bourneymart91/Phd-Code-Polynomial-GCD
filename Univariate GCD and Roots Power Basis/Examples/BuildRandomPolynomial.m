function [root_mult_arr] = BuildRandomPolynomial(m,int_low,int_up)
% Get the a array of roots and their corresponding multiplicities for a
% random polynomial, where the roots are given within an interval [a,b]
%
% Inputs.
%
% m : Degree of random polynomial
%
% int_low : Lower bound for roots interval
%
% int_up : Upper bound for root interval
%
% Example: BuildRandomPolynomial(10,-5,5)

global SETTINGS

% Set the upper and lower bound
a = int_low;
b = int_up;

% Get a multiplicity structure for polynomial f
    
    % Generate a probability array, so that we produce some roots of high
    % multiplicity.
    prob_arr = zeros(1,m);
    for i = 1:1:m
        prob_arr(i) = i./ nchoosek(m+1,2);
    end
    
    prob_arr = fliplr(prob_arr);
    
    rng(SETTINGS.SEED);
    
    % Get the multiplicity structure of d.
    total = 0;
    i = 1;
    while total < m
        r = rand;
        prob = prob_arr;
        x = sum(r >= cumsum([0, prob]));
        if (total + x) <= m
            mult_arr_d(i) = x;
            total = total + x;
            i = i+1;
        end
    end
    
    % get the number of roots of d
    num_roots_t = length(mult_arr_d);
    
    
    % Get the roots
    % Get a set of unique roots
    % the 1000 and 1000 contain the roots to the unit interval
    detail = 100;
    format 'long';


    
    roots = a + randperm(detail,num_roots_t)./(detail./(b-a));
    
    root_mult_arr = [roots' mult_arr_d'];
end
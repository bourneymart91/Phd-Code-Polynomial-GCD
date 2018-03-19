function [] = test3

% Compute b_{i,j} given i,j, and n

n = 3;

i = 2;
j = 1;

syms a_00 a_01 a_02 a_03
syms a_10 a_11 a_12 a_13
syms a_20 a_21 a_22 a_23
syms a_30 a_31 a_32 a_33

% b_{i,j} = 

mat_coef_f = ...
    [
         a_00 a_01 a_02 a_03
         a_10 a_11 a_12 a_13
         a_20 a_21 a_22 a_23
         a_30 a_31 a_32 a_33
    ];

bij = 0;

for t = 0:1:j
    for s = 0:1:i
        bij = bij + ...
            (nchoosek(i,s) * nchoosek(j,t) / Trinomial(n,s,t)) ...
            *mat_coef_f(s+1,t+1);
    end
end

display(i)
display(j)
display(n)
display(bij)

end


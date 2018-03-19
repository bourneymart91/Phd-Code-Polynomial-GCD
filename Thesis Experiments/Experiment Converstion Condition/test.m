
% Get a set of roots and multiplicities
mult_root_array = ...
    [
    1.0687 2;
    1.2547 3; 
    1.2783 4;
    2.9871 5; 
    ];




f_bb = GetCoefficients_Bernstein(mult_root_array);
f_pwr = GetCoefficients_Power(mult_root_array);


bb_rep = PowerToBernstein(f_pwr)
pwr_rep = BernsteinToPower(f_bb)


% Get difference between f_bb and bb_rep
norm(f_bb - bb_rep)


% Get difference between f_pwr and pwr_rep
norm(f_pwr - pwr_rep)



function fx = Examples_Roots_FromCoefficients(ex_num)

addpath(genpath('../Examples'));

% Get the symbolic factors of f(x) and corresponding multiplicities
f_root_sym_mult_array = Roots_Examples_Univariate(ex_num);

display(f_root_sym_mult_array)

% Get the coefficients of f(x) as a vector.
fx = GetCoefficientsFromSymbolicRoots(f_root_sym_mult_array);


end


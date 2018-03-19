

arr_noise = {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12};

for i = 1:1:length(arr_noise)

    try
        o_gcd_Univariate_3Polys('5', arr_noise{i}, 1e-12, 'Geometric Mean Matlab Method', true, 'None', 'None')
    catch
    end
end


for i = 1:1:length(arr_noise)

    try
        o_gcd_Univariate_3Polys('5', arr_noise{i}, 1e-12, 'Geometric Mean Matlab Method', false, 'None', 'None')
    catch
    end
end


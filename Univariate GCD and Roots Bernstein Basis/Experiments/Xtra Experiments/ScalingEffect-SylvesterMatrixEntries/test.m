function [] = test(m, n_k)
% test(m, n_k)
%
% % Inputs
%
% m : Degree of polynomial f(x)
%
% n_k : Degree of polynomial v(x)
%
% Example
%
% >> test(5, 7)


% Initialise vector to store sum of scalars for each a_{i}
v_scalar_sum = zeros(m+1,1);
v_scalar_sum2_numerator = zeros(m+1,1);
v_scalar_sum2_numerator_part1 = zeros(m+1,1);
v_scalar_sum2_numerator_part2 = zeros(m+1,1);

for i = 0:1:m
    temp_sum = 0;
    temp_sum2_numerator = 0;
    temp_sum2_numerator_part1 = 0;
    temp_sum2_numerator_part2 = 0;
    
    for j = 0:1:n_k
        temp_sum = temp_sum + (nchoosek(m,i)  * nchoosek(n_k,j)  ./ nchoosek(m+n_k,i+j));
        
        temp_sum2_numerator = temp_sum2_numerator + (nchoosek(i+j,i)*nchoosek(m+n_k-i-j,m-i));
        
        temp_sum2_numerator_part1 = temp_sum2_numerator_part1 + (nchoosek(i+j,i));
        temp_sum2_numerator_part2 = temp_sum2_numerator_part2 + (nchoosek(m+n_k-i-j,m-i));
        
    end
    
    
    v_scalar_sum(i+1) = temp_sum;
    v_scalar_sum2_numerator(i+1) = temp_sum2_numerator;
    
    v_scalar_sum2_numerator_part1(i+1) = temp_sum2_numerator_part1;
    v_scalar_sum2_numerator_part2(i+1) = temp_sum2_numerator_part2;
    
end

display(v_scalar_sum)
display(v_scalar_sum2_numerator)

display(v_scalar_sum2_numerator_part1)
display(v_scalar_sum2_numerator_part2)

denominator = nchoosek(m+n_k,m);
display(denominator)


display(v_scalar_sum2_numerator(1)./nchoosek(m+n_k,m))

end
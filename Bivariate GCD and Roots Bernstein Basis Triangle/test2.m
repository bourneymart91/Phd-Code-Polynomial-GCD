

m = 5;
n_k = 7;


first_binomial = zeros(n_k+1, n_k+1);
second_binomial = zeros(n_k+1, n_k+1);
third_binomial = zeros(n_k+1, n_k+1);
fourth_binomial = zeros(n_k+1, n_k+1);

prod_first_binomial = zeros(m+1,m+1);
prod_second_binomial = zeros(m+1,m+1);
prod_third_binomial = zeros(m+1,m+1);

for k1 = 0:1:m
    for i1 = k1:-1:0
        
        i2 = k1 - i1;
        
        first_binomial = zeros(n_k+1, n_k+1);
        second_binomial = zeros(n_k+1, n_k+1);
        third_binomial = zeros(n_k+1, n_k+1);
        fourth_binomial = zeros(n_k+1, n_k+1);
        
        coeff_multiplier = zeros(n_k+1, n_k+1);
        
        for k2 = 0:1:n_k
            for j1 = k2:-1:0
                
                j2 = k2 - j1;
                
                first_binomial(j1+1, j2+1) = nchoosek(i1 + j1, i1);
                second_binomial(j1+1, j2+1) = nchoosek(i2 + j2, i2);
                third_binomial(j1+1,j2+1) = nchoosek(m+n_k - i1 - i2 - j1 - j2, m - i1 - i2);
                
                fourth_binomial(j1+1,j2+1) = nchoosek(m+n_k, n_k);
                %fourth_binomial(j1+1,j2+1) = 1;
                
                coeff_multiplier(j1+1,j2+1) = first_binomial(j1+1,j2+1) * second_binomial(j1+1,j2+1) * third_binomial(j1+1,j2+1) ; %/ fourth_binomial(j1+1,j2+1);
                
                
                
            end
        end
        
        display(coeff_multiplier)
        
        prod_first_binomial(i1+1, i2+1) = prod(first_binomial(first_binomial~=0));
        prod_second_binomial(i1+1, i2+1) = prod(second_binomial(second_binomial~=0));
        prod_third_binomial(i1+1, i2+1) = prod(third_binomial(third_binomial~=0));
        
        myprod(i1+1, i2+1) = prod(coeff_multiplier(coeff_multiplier~=0));
        
    end
end

display(prod_first_binomial)
display(prod_second_binomial)
display(prod_third_binomial)


myprod2 = prod(myprod(myprod~=0))

display(first_binomial)
display(second_binomial)
display(third_binomial)
display(fourth_binomial)

display(myprod)

function [] = GeometricMeanCols(m,n_k)

for j = 0:1:n_k
    
    temp_prod = 1;
    for i = 0:1:m
        
        vColumnValues(i+1) = nchoosek(m,i) * nchoosek(n_k,j) ./ nchoosek(m+n_k,i+j);
        
    end
    
    vColumnMean(j+1) = geomean(vColumnValues);
    
end

figure()
hold on
plot(vColumnMean)
hold off

1./vColumnMean



end
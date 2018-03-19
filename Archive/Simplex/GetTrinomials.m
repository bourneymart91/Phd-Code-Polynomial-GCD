function Tri_f = GetTrinomials(m)


% Build matrix of trinomials
mat = zeros(m+1,m+1)

for i = 0:1:m
    for j = 0:1:m
        if i+j <= m
            mat(i+1,j+1) = Trinomial(m,i,j);
        end
    end
end
    
Tri_f_vec = GetAsVector(mat);

% Remove zeros
nCoefficients = nchoosek(m+2,2);

Tri_f = Tri_f_vec(1:nCoefficients);

end
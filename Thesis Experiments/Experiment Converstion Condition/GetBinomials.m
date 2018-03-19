function [Bi_m] = GetBinomials(m)

Bi_m = zeros(m+1,1);

for i = 0:1:m
    Bi_m(i+1) = nchoosek(m,i);
end

end
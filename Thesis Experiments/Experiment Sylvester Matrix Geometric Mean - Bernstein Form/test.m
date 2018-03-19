
m = 6

n = 4
k = 2



prod1 = 1;
prod2 = 1;

for i = 0:1:m
    prod1 = prod1 .* nchoosek(n-k+1+i,i)
    prod2 = prod2 .* nchoosek(m+n-k+1-i,m-i)
end

display(prod1)
display(prod2)

prod1 .^ (1./((n-k+1).*(m+1)))
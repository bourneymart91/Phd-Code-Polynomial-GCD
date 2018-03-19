function nck = Comb1(n,k)

if (k == 0) || (k==n)
    nck = 1; 
else
    nck = Comb1(n-1,k-1) + Comb1(n-1,k);
end

end
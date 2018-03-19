function nCk = lnnchoosek(n,k)

nCk = sum(log10(1:1:n)) ...
    - sum(log10(1:1:k)) ...
    - sum(log10(1:1:(n-k)));

end
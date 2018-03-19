function nCk = Trinomial(m,i,j)

nCk = factorial(m)./(factorial(i)*factorial(j)*factorial(m-i-j));

end
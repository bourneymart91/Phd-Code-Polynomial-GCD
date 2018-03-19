function [] = Rearranged_test1(m,n_k)
% Testing the relationship between product of 
%
% 1. Product of a_{i} in S_{k} and product of a_{i} in S_{k-1}
%
% 2. Product of \binom{i+j,j} in S_{k} and in S_{k-1}
%
% 3. Product of \binom{m+n-k-i-j}{m-i} in S_{k} and S_{k-1}
%
% 4. Product of \frac{1} {\binom{m+n-k}{m}} in S_{k} and S_{k-1}

% Part One : Coefficients

% Part Two : product of nchoosek(i+j,j)

% %
% 
% %
% %
% Test Part Two
i = 1;
test1a = 1;
for j = 0:1:n_k
   test1a = test1a .* nchoosek(i+j,j);
end

test1b = 1;
for j = 0:1:n_k+1
    test1b = test1b .* nchoosek(i+j,j);
end

test_Part_Two_equivalenct = (test1b * 1./(nchoosek(i+n_k+1,i)) ) - test1a;
% % 
% % 
% Test Part Three

test3a = 1;
for j = 0:1:n_k
    test3a = test3a * nchoosek(m+n_k-i-j,m-i);
end

test3b = 1;
for j = 0:1:n_k+1
    test3b = test3b * nchoosek(m+n_k-i-j+1,m-i);
end

test_Part_Three_Equivalence = (test3b * (1./ nchoosek(m+n_k-i+1,m-i)) ) - test3a;

% %
% %
% % 
% Test Part Four

test4a =1;
for j = 0:1:n_k
    test4a= test4a .* (1./ nchoosek(m+n_k,m));
end

test4b = 1;
for j = 0:1:n_k+1
    test4b = test4b .* ...
        ( 1./ nchoosek(m+n_k+1,m) );    
end

test_Part4_Equivalence =  (test4b * nchoosek(m+n_k,m) ./ ...
    (((n_k+1)./(m+n_k+1))^(n_k+2))) - test4a;




end
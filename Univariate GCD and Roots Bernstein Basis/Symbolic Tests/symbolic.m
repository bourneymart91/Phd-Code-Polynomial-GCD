function A = symbolic(n)

% A = sym(zeros(n, n));
% for i = 1:n
%     for j = 1:n
%         A(i, j) = sym(sprintf('a%d%d', i, j));
%     end
% end



f = [0 1 2 3 4 5 6];
g = [0 1 2 3];

% get degree of polynomial f.
m = length(f) -1;
% get degree of polynomial g.
n = length(g) -1;
k = 1;

A = sym(zeros(m+n-k+1, n-k+1));
for j = 0:1:n-k
    for i = j:1:j+m
        A(i+1,j+1) = sym(sprintf('a%d', i-j)) .* nchoosek(m,i-j);
    end
end
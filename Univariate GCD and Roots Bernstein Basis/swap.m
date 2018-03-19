function [a,b] = swap(a,b)

temp = b;
% Replace b with a value
b = a;
% Replace a with temp value
a = temp;

end
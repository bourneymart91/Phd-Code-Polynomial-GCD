function Y = BuildY_STLN(x,m,n,k)
% Build the matrix Y_{k}, where Y_{k}z = E_{t}(z)x, 
% Similarly Y_{t}[f;g] = S_{t}(f,g) * x
%
% Used in STLN function
%
% % Inputs.
%
% x : The solution to S_{t}(f,g) x = c_{t}
% 
% m : Degree of polynomial f(x)
% 
% n : Degree of polynomial g(x)
%
% k : Index of kth Sylvester Subresultant matrix
%
% % Outputs.
%
% Y : Matrix such that Y_{k} such that Y_{k}(x1,x2)z = S(z1,z2) x

% Split the vector x
x1 = x(1:n-k+1);
x2 = x(n-k+2:end);

% Build the matrices Y1 and Y2 - the two partitions of Y.
Y1 = BuildT1(x1,m);
Y2 = BuildT1(x2,n);

% Build Y
Y = [Y1 Y2];


end
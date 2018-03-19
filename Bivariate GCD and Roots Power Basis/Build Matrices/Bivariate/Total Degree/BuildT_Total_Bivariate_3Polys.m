function Sk = BuildT_Total_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, k)
% Build the kth Sylvester matrix where polynomials f(x,y) and g(x,y) are
% given in terms of their total degree.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of polynomial h(x,y)
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% o : (Int) Total degree of h(x,y)
%
% k : Total degree of d(x,y) and index of the kth Sylvester subresultant
% matrix to be built.
%
% % Outputs.
%
% Sk : (Matrix) The kth Sylvester matrix S_{k}

% Build the partitions of the Sylvester matrix
T1 = BuildT1_Total_Bivariate(fxy, m, n-k);
T2 = BuildT1_Total_Bivariate(fxy, m, o-k);
T3 = BuildT1_Total_Bivariate(gxy, n, m-k);
T4 = BuildT1_Total_Bivariate(hxy, o, m-k);

% Build the kth Sylvester subresultant matrix.
diagonal = blkdiag(T1, T2);
column = [T3 ; T4];
Sk = [diagonal column]; 

end
function Sk = BuildT_Total_Bivariate_2Polys(fxy,gxy,m,n,k)
% Build the kth Sylvester matrix where polynomials f(x,y) and g(x,y) are
% given in terms of their total degree.
%
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% k : (Int) Total degree of d(x,y) and index of the kth Sylvester subresultant
% matrix to be built.
%
% % Outputs.
%
% Sk : (Matrix) The kth Sylvester matrix S_{k}

% Build the partition T_{n-k}(f)
Tf = BuildT1_Total_Bivariate(fxy, m, n-k);

% Build the partiton T_{m-k}(g)
Tg = BuildT1_Total_Bivariate(gxy, n, m-k);

% Build the kth Sylvester subresultant matrix.
Sk = [Tf Tg]; 

end
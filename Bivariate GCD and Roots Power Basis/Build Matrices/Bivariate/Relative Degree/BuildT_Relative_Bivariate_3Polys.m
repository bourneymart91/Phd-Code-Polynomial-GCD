function Sk = BuildT_Relative_Bivariate_3Polys(fxy, gxy, hxy, k1,k2)
% Given two input polynomials f(x,y) and g(x,y), build the (k1,k2)-th
% Sylvester subresultant.
%
% Inputs
%
% fxy : (Matrix) Coefficients of input polynomial f(x,y)
% 
% gxy : (Matrix) Coefficients of input polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of input polynomial h(x,y)
%
% k1 : (Int) Degree of common divisor with respect to x.
%
% k2 : (Int) Degree of common divisor with respect to y.
%
% Outputs.
%
% Sk : (Matrix) The Sylvester Subresultant S_{k_{1},k_{2}}


% Get degrees of polynomial f(x,y), g(x,y) and h(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);

% Build the partitions of the Sylvester subresultant matrix S_{k1,k2}(f,g)
T1 = BuildT1_Relative_Bivariate_Version1(fxy, n1-k1, n2-k2);
T2 = BuildT1_Relative_Bivariate_Version1(fxy, o1-k1, o2-k2);
T3 = BuildT1_Relative_Bivariate_Version1(gxy, m1-k1, m2-k2);
T4 = BuildT1_Relative_Bivariate_Version1(hxy, m1-k1, m2-k2);

% Build the Sylvester subresultant matrix
diagonal = blkdiag(T1,T2);
column = [T3; T4];
Sk = [diagonal column];


end
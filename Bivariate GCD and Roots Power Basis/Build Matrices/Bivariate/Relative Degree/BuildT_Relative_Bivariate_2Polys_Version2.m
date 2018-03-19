function T = BuildT_Relative_Bivariate_2Polys_Version2(fxy, gxy, k1, k2)
% Build the Sylvester matrix 
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
% 
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% k1 : (Int) Degree of common divisor d(x,y) with respect to x
%
% k2 : (Int) Degree of common divisor d(x,y) with respect to y
%
%
% % Outputs
%
% T : Bivariate Sylvester matrix [T(f) T(g)]


% Get degree of f(x,y) and g(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

% Build the partitions of the Sylvester matrix
T1 = BuildT1_Relative_Bivariate_Version2(fxy, n1-k1, n2-k2);
T2 = BuildT1_Relative_Bivariate_Version2(gxy, m1-k1, m2-k2);

T = [T1 T2];


    
end

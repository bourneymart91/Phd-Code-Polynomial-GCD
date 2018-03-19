function T = BuildT_3Polys_3Eqns(fxy, gxy, hxy, m, n, o, k)
%
% % Inputs
%
% fxy : (Matrix)
%
% gxy : (Matrix)
%
% hxy : (Matrix)
%
% m : (Int)
%
% n : (Int)
%
% o : (Int)
% 
% k : (Int)


% Build the matrix T
T1 = BuildT1(fxy, m, n - k);
T3 = BuildT1(gxy, n, m - k);
T5 = BuildT1(fxy, m, o - k);
T6 = BuildT1(hxy, o, m - k);
T7 = BuildT1(hxy, o, n - k);
T8 = BuildT1(gxy, n, o - k);


nRows = nchoosek(m + n - k + 2, 2);
nCols = nchoosek(o - k + 2, 2);
T2 = zeros(nRows, nCols);

nRows = nchoosek(m + o - k + 2, 2);
nCols = nchoosek(n - k + 2, 2);
T4 = zeros(nRows, nCols);

nRows = nchoosek(n + o - k + 2, 2);
nCols = nchoosek(m - k + 2, 2);
T9 = zeros(nRows, nCols);


T = [T1 T2 T3; T4 T5 T6; T7 -T8 T9];


end
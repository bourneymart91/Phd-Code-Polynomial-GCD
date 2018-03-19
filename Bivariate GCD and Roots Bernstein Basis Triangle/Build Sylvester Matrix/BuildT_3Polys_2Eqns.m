function T = BuildT_3Polys_2Eqns(fxy, gxy, hxy, m, n, o, k)
% 
% % Build the matrix T
% T1_f = BuildT1(fxy, m, n - k);
% T2_f = BuildT1(fxy, m, o - k);
% T3_g = BuildT1(gxy, n, m - k);
% T4_h = BuildT1(hxy, o, m - k);


T1 = BuildT1(fxy, m, n - k);
T2 = zeros(nchoosek(m + n - k + 2, 2), nchoosek(o - k + 2, 2));
T3 = BuildT1(gxy, n, m - k);

T4 = zeros(nchoosek(m + o - k + 2, 2), nchoosek(n - k + 2, 2));
T5 = BuildT1(fxy, m, o - k);
T6 = BuildT1(hxy, o, m - k);

T = [T1 T2 T3 ; T4 T5 T6];


end
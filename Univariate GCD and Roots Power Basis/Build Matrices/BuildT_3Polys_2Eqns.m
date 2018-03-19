function T = BuildT_3Polys_2Eqns(fx, gx, hx, k)
% Build the Sylvester Subresultant matrix S_{k}(f,g)
%
% % Inputs
%
% fx : Coefficients of the polynomial f(x)
%
% gx : Coefficients of the polynomial g(x)
%
% hx : Coefficients of the polynomial h(x)
%
% k : Index of Sylvester Subresultant matrix to be constructed.


% Get degree of polynomial f(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

T1 = BuildT1(fx, n - k);
T2 = zeros(m + n - k + 1, o - k  + 1);
T3 = BuildT1(gx, m - k);
T4 = zeros(m + o - k + 1, n - k + 1);
T5 = BuildT1(fx, o - k);
T6 = BuildT1(hx, m - k);

T = [T1 T2 T3; T4 T5 T6];


block_row_1 = [T1 T2 T3];
block_row_2 = [T4 T5 T6];

% Joab's Ratio

v_br1 = block_row_1(block_row_1 ~= 0);
v_br2 = block_row_2(block_row_2 ~= 0);


max_br1 = max(max(abs(v_br1)));
min_br1 = min(min(abs(v_br1)));

ratio_br1 = max_br1 ./ min_br1;

max_br2 = max(max(abs(v_br2)));
min_br2 = min(min(abs(v_br2)));

ratio_br2 = max_br2 ./ min_br2;

LineBreakSmall()
fprintf('r1 : %e \n', ratio_br1);
fprintf('r2 : %e \n', ratio_br2);
fprintf('Ratio : Block row 1 : Block row 2 = %e \n', ratio_br1 ./ ratio_br2);
LineBreakSmall()


% % Martin's Ratio
vT1 = T1(T1~=0);
vT3 = T3(T3~=0);
vT5 = T5(T5~=0);
vT6 = T6(T6~=0);


max_T1 = max(max(abs(vT1)));
max_T3 = max(max(abs(vT3)));
max_T5 = max(max(abs(vT5)));
max_T6 = max(max(abs(vT6)));

min_T1 = min(min(abs(vT1)));
min_T3 = min(min(abs(vT3)));
min_T5 = min(min(abs(vT5)));
min_T6 = min(min(abs(vT6)));



ratio_T1 = max_T1 ./ min_T1;
ratio_T3 = max_T3 ./ min_T3;
ratio_T5 = max_T5 ./ min_T5;
ratio_T6 = max_T6 ./ min_T6;

LineBreakLarge()
fprintf('')
fprintf('Max / Min T1 : %e \n', ratio_T1)
fprintf('Max / Min T3 : %e \n', ratio_T3)
fprintf('Max / Min T5 : %e \n', ratio_T5)
fprintf('Max / Min T6 : %e \n', ratio_T6)
LineBreakLarge()
fprintf('')
fprintf('Ratio Columns in First Row-Partition : %e \n', ratio_T1 ./ ratio_T3);
fprintf('Ratio Columns in Second Row-Partition : %e \n', ratio_T5 ./ ratio_T6);
fprintf('Ratio Rows in Third Column-Partition : %e \n', ratio_T3 ./ ratio_T6);
LineBreakLarge()







end
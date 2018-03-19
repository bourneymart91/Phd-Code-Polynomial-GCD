function [] = analyseResults_Roots(ex_num)
%
% % Inputs
%
% ex_num : (String)

T = readtable("Results_o_Roots.dat");

% Filter based on example number
filteredTable = T(T.EX_NUM == ex_num, : )

end
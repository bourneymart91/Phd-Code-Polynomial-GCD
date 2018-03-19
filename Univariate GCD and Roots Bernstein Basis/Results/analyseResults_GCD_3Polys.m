function [] = analyseResults_GCD_3Polys(ex_num)
%
% % Inputs
%
% ex_num : (String) Example Number


% Get table of results
T = readtable("Results_o_GCD_3Polys.dat");

% Filter based on example number
filteredTable = T(strcmp(T.EX_NUM, ex_num), : )

end
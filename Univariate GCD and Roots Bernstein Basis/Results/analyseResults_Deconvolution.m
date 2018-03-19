function [] = analyseResults_Deconvolution(ex_num)
%
% % Inputs
%
% ex_num : (String) Example Number


% Get table of results
T = readtable("Results_o_deconvolutions.dat");

% Filter based on example number
filteredTable = T(T.EX_NUM == ex_num, : )

end
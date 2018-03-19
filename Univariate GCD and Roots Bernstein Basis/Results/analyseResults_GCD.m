function [] = analyseResults_GCD(ex_num)





T = readtable("Results_o_gcd.dat");




filteredTable = T(strcmp(T.EX_NUM , num2str(ex_num)), : )



end
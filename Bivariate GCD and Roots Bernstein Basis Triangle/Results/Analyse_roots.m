function [] = Analyse_roots(ex_num)


% 19-Jan-2018 13:58:05,1,1.489497e-10,2.106368e-10,1.671762e-10,Geometric Mean Matlab Method,1,None,0,None,0,1.000000e-10,1.000000e-08,DTQ,Minimum Singular Values 
%
%     DATE,                 (String)
%     EX_NUM,               (String)
%     ERROR_FXY,            (Float)
%     ERROR_HXY,            (Float)
%     ERROR_WXY,            (Float)
%     MEAN_METHOD,          (String)
%     BOOL_ALPHA_THETA,     (String)
%     LOW_RANK_APPROX_METHOD, (String)
%     LRA_ITE,                  (Int)
%     APF_METHOD,               (String)
%     APF_ITE,                  (Int)
%     error_min,                (Float)
%     error_max,                (Float)
%     SUBRESULTANT_FORMAT,      (String)
%     RANK_REVEALING_METRIC     (String)
%     DECONVOLUTION_METHOD_HXY (String)
%     DECONVOLUTION_METHOD_WXY (String)




T = readtable("Results_o_roots.dat");



% Filter based on example number
filteredTable = T(T.EX_NUM == ex_num, : )

end
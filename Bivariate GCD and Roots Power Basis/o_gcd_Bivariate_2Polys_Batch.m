function [] = o_gcd_Bivariate_2Polys_Batch
% Performs a batch of GCD computations for a variety of settings

arrExNum = {'1','2','3','4','5','6','7','8','9','10','11'};
%ex_num_arr = {'1','2'};
arrEmin = {1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12};
%emin_arr = {1e-8,1e-10};
arrEmax = {1e-12};
arrMean_method = {'Geometric Mean Matlab Method', 'None'};
arrBool_alpha_theta = {true, false};
arrLow_rank_approx_method = {'None', 'Standard STLN', 'Standard SNTLN'};
arrAPFMethod = {'None'};
arrDegreeMethod = {'Relative', 'Total', 'Both'};
arrRankRevealingMetric = {'Minimum Singular Values'} ; %'R1 Row Norms' not complete

parfor i1 = 1:1:length(arrExNum)
    
    for i2 = 1:1:length(arrEmin)
        
        for i3 = 1:1:length(arrMean_method)
            
            for i4 = 1:1:length(arrBool_alpha_theta)
                
                for i5 = 1:1:length(arrLow_rank_approx_method)
                    
                    for i6 = 1:1:length(arrAPFMethod)
                        
                        for i7 = 1:1:length(arrRankRevealingMetric)
                            
                            
                            ex_num = arrExNum{i1};
                            emin = arrEmin{i2};
                            emax = arrEmax{1};
                            mean_method = arrMean_method{i3};
                            bool_alpha_theta = arrBool_alpha_theta{i4};
                            low_rank_approx_method = arrLow_rank_approx_method{i5};
                            apf_method = arrAPFMethod{i6};
                            rank_revealing_metric = arrRankRevealingMetric{i7};
                            
                            
                            for i8 = 1:1:length(arrDegreeMethod)
                                
                                degree_method = arrDegreeMethod{i8};
                                
                                
                                try
                                    
                                    close all;
                                    clc;
                                    o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, degree_method, rank_revealing_metric)
                                    fileId = fopen('log.txt','a')
                                    fprintf(fileId,'%s','success \n');
                                    fclose(fileId);
                                    
                                catch err
                                    
                                    fileId = fopen('log.txt','a')
                                    fprintf(fileId,'%s \n\n\n',getReport(err));
                                    fclose(fileId);
                                    
                                end
                                
                            end
                        end
                    end
                end
            end
        end
    end
end

end
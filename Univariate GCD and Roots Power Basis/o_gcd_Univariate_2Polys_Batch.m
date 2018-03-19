function [] = o_gcd_Univariate_2Polys_Batch()
% Perform a batch of gcd examples
%
%
% % Example
%
% o_gcd_batch() 

arrExampleNumber = {'1','2','3','4','5','6','7','8','9','10'};
arrEmin = {1e-8, 1e-10, 1e-12};
arrEmax = {1e-12};
arrMeanMethod = {'Geometric Mean Matlab Method', 'None'};
arrBoolAlphaTheta = {true, false};
arrLowRankApproxMethod = {'Standard STLN', 'Standard SNTLN', 'None'};
%apf_method_arr = {'None','Standard APF Nonlinear','Standard APF Linear'};
arrAPFMethod = {'None'};


parfor i1 = 1:1:length(arrExampleNumber)
    
    ex_num = arrExampleNumber{i1};
    
    for i2 = 1:1:length(arrEmin)
        
        emin = arrEmin{i2};
        
        for i3 = 1:1:length(arrEmax)
            
            emax = arrEmax{i3};
            
            for i4 = 1:1:length(arrMeanMethod)
                
                mean_method = arrMeanMethod{i4};
                
                for i5 = 1:1:length(arrBoolAlphaTheta)
                    
                    bool_alpha_theta = arrBoolAlphaTheta{i5};
                    
                    for i6 = 1:1:length(arrLowRankApproxMethod)
                        
                        low_rank_approx_method = arrLowRankApproxMethod{i6};
                        
                        for i7 = 1:1:length(arrAPFMethod)
                            
                            apf_method = arrAPFMethod{i7};
                            
                            try
                                close all;
                                clc;
                                o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method)
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

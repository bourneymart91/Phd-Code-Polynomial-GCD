function [] = o_gcd_Bivariate_3Polys_batch
% Performs a batch of GCD computations for a variety of settings

ex_num_arr = {'1','2','3','4','5','6','7','8','9','10','11'};
%ex_num_arr = {'1','2'};
emin_arr = {1e-7,1e-8,1e-9,1e-10,1e-11,1e-12};
emin_arr = {1e-8,1e-10};
emax_arr = {1e-12};
mean_method_arr = {'Geometric Mean Matlab Method','None'};
bool_alpha_theta_arr = {'y','n'};
low_rank_approx_method_arr = {'None','Standard STLN','Standard SNTLN'};
apf_method_arr = {'None'};
degree_method_arr = {'Relative','Total','Both'};

for i1 = 1:1:length(ex_num_arr)
    
    for i2 = 1:1:length(emin_arr)
        
        for i3 = 1:1:length(mean_method_arr)
            
            for i4 = 1:1:length(bool_alpha_theta_arr)
                
                for i5 = 1:1:length(low_rank_approx_method_arr)
                    for i6 = 1:1:length(apf_method_arr)
                        
                        ex_num = ex_num_arr{i1};
                        emin = emin_arr{i2};
                        emax = emax_arr{1};
                        mean_method = mean_method_arr{i3};
                        bool_alpha_theta = bool_alpha_theta_arr{i4};
                        low_rank_approx_method = low_rank_approx_method_arr{i5};
                        apf_method = apf_method_arr{i6};
                        
                        
                        parfor i7 = 1:1:length(degree_method_arr)
                            
                            degree_method = degree_method_arr{i7};
                            
                            
                            try
                                close all;
                                clc;
                                o_gcd_Bivariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, degree_method)
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
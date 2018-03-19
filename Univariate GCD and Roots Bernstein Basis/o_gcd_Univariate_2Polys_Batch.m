function [] = o_gcd_Univariate_2Polys_Batch()
% Perfome a batch of univariate gcd computations with varying input
% parameters.

arr_ex_num = {'1','2','3','4','5','6','7','8','9','10'};
arr_bool_alpha_theta = {true, false};
arr_emin = {1e-08, 1e-09, 1e-10, 1e-11, 1e-12};
arr_low_rank_approx_method = {'Standard SNTLN','Standard STLN','None'};
arr_apf_method = {'None'};
arr_mean_method_arr = {'Geometric Mean Matlab Method', 'None', 'Arithmetic Mean'};
arr_Sylvester_matrix_variant = {'T', 'DT', 'DTQ', 'TQ', 'DTQ Denominator Removed'};
arr_rank_revealing_metric = {'Minimum Singular Values', 'Max/Min Singular Values', 'R1 Row Norms', 'R1 Row Diagonals', 'Residuals'};


arr_bool_log = {false};
arr_gcd_coefficient_method = {'ux and vx'};

global SETTINGS



for i1 = 1:1:length(arr_bool_log)
    
    SETTINGS.BOOL_LOG = arr_bool_log{i1};
    
    for i2 = 1:1:length(arr_gcd_coefficient_method)
        
        SETTINGS.GCD_COEFFICIENT_METHOD = arr_gcd_coefficient_method{i2};
        
        % Changing example number
        parfor i3 = 1 : 1 : length(arr_ex_num)
            for i4 = 1 : 1 : length(arr_emin)
                for i5 = 1 : 1 : length(arr_low_rank_approx_method)
                    for i6 = 1 : 1 : length(arr_bool_alpha_theta)
                        for i7 = 1 : 1 : length(arr_mean_method_arr)
                            for i8 = 1 : 1 : length(arr_Sylvester_matrix_variant)
                                for i9 = 1 : 1 : length(arr_apf_method)
                                    for i10 = 1 : 1 : length(arr_rank_revealing_metric)
                                        
                                        ex_num = arr_ex_num{i3};
                                        emin = arr_emin{i4};
                                        emax = 1e-12;
                                        low_rank_approx_method = arr_low_rank_approx_method{i5};
                                        bool_alpha_theta = arr_bool_alpha_theta{i6};
                                        mean_method = arr_mean_method_arr{i7};
                                        sylvester_matrix_variant = arr_Sylvester_matrix_variant{i8};
                                        apf_method = arr_apf_method{i9};
                                        rank_revealing_metric = arr_rank_revealing_metric{i10};
                                        
                                        
                                        
                                        file_name = 'log_GCD_computations.txt'
                                        try
                                            close all;
                                            clc;
                                            o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_matrix_variant, rank_revealing_metric)
                                            fileId = fopen(file_name,'a')
                                            fprintf(fileId,'%s %s \n', datetime('now'), 'Success' );
                                            fclose(fileId);
                                            
                                        catch err
                                            
                                            fileId = fopen(file_name,'a')
                                            fprintf(fileId,'%s %s \n\n\n', datetime('now') , getReport(err));
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
    
    
end







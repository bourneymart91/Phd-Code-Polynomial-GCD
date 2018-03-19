function [] = o_roots_Bivariate_batch()
% o_roots_Bivariate_batch()
%
%

ex_num_arr = {'1', '2', '3', '4'};
emin = 1e-12;
emax_arr = {1e-8, 1e-10, 1e-12};
mean_method_arr = {'Geometric Mean Matlab Method', 'None'};
bool_alpha_theta_arr = {true, false};
low_rank_approx_method_arr = {'None', 'Standard STLN', 'Standard SNTLN'};
degree_method_arr = {'Total', 'Relative', 'Both'};
arrRankRevealingMetric = {'Minimum Singular Values'};


parfor i1 = 1:1:length(ex_num_arr)
    for i2 = 1:1:length(emax_arr)
        for i3 = 1:1:length(mean_method_arr)
            for i4 = 1:1:length(bool_alpha_theta_arr)
                for i5 = 1:1:length(low_rank_approx_method_arr)
                    for i6 = 1:1:length(degree_method_arr)
                        for i7 = 1:1:length(arrRankRevealingMetric)
                            for i8 = 1:1:length(arrDeconvolutionMethod_hxy)
                                for i9 = 1:1:length(arrDeconvolutionMethod_wxy)
                            
                            
                                    ex_num = ex_num_arr{i1};
                                    emax = emax_arr{i2};
                                    mean_method = mean_method_arr{i3};
                                    bool_alpha_theta = bool_alpha_theta_arr{i4};
                                    low_rank_approx_method = low_rank_approx_method_arr{i5};
                                    degree_method = degree_method_arr{i6};
                                    rank_revealing_metric = arrRankRevealingMetric{i7};
                                    deconvolution_method_hxy = arrDeconvolutionMethod_hxy{i8}
                                    deconvolution_method_wxy = arrDeconvolutionMethod_wxy{i9}
                                    
                                    apf_method = 'None';

                                    o_roots_Bivariate(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, degree_method, rank_revealing_metric, deconvolution_method_hxy, deconvolution_method_wxy);
                                
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
function [] = o_roots_Univariate_batch



arrExampleNumber = {'1','2','3','4','5','6','7','8','9','10'};
arrEmax = {1e-10, 1e-8, 1e-6};

arrMeanMethod = {'None', 'Geometric Mean Matlab Method'};

arrBoolAlphaTheta = {true, false};
arrLowRankApproxMethod = {'None', 'Standard STLN', 'Standard SNTLN'};
arrAPFMethod = {'None'};
arrSylvesterMatrixVariant = {'DTQ'};%, 'TQ', 'DT', 'T'};
arrDeconvolutionMethod_hx = {'Separate', 'Batch', 'Batch Constrained', 'Batch Constrained With STLN', 'Batch With STLN'};
arrDeconvolutionMethod_wx = {'Separate', 'Batch'};
arrRankRevealingMetric = {'Minimum Singular Values'};%, 'Max/Min Singular Values', 'R1 Row Norms', 'R1 Row Diagonals', 'Residuals'};
arrDeconvolutionPreproc = {false, true};


parfor i1 = 1:1:length(arrExampleNumber)
    for i2 = 1:1:length(arrEmax)
        for i3 = 1:1:length(arrMeanMethod)
            for i4 = 1:1:length(arrBoolAlphaTheta)
                for i5 = 1:1:length(arrLowRankApproxMethod)
                    for i6 = 1:1:length(arrAPFMethod)
                        for i7 = 1:1:length(arrSylvesterMatrixVariant)
                            for i8 = 1 : 1 : length(arrDeconvolutionMethod_hx)
                                for i9 = 1:1:length(arrDeconvolutionMethod_wx)
                                    for i10 = 1:1:length(arrRankRevealingMetric)
                                        for i11 = 1:1:length(arrDeconvolutionPreproc)
                                            
                                            ex_num = arrExampleNumber{i1};
                                            emin = 1e-12;
                                            
                                            emax = arrEmax{i2};
                                            mean_method = arrMeanMethod{i3};
                                            bool_alpha_theta = arrBoolAlphaTheta{i4};
                                            low_rank_approx_method = arrLowRankApproxMethod{i5};
                                            apf_method = arrAPFMethod{i6};
                                            sylvester_matrix_variant = arrSylvesterMatrixVariant{i7};
                                            deconvolution_method_hx =  arrDeconvolutionMethod_hx{i8};
                                            deconvolution_method_wx = arrDeconvolutionMethod_wx{i9};
                                            rank_revealing_metric = arrRankRevealingMetric{i10};
                                            deconvolution_preproc = arrDeconvolutionPreproc{i11};
                                            
                                            file_name = 'log_roots.txt';
                                            try
                                                close all;
                                                clc;
                                                o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
                                                    low_rank_approx_method, apf_method, sylvester_matrix_variant, rank_revealing_metric, deconvolution_method_hx, deconvolution_method_wx, deconvolution_preproc)
                                                fileId = fopen(file_name,'a');
                                                fprintf(fileId,'%s %s \n', datetime('now'), 'Success' );
                                                fclose(fileId);
                                                
                                            catch err
                                                
                                                fileId = fopen(file_name,'a');
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
    
    
end



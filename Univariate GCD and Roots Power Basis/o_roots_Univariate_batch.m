

function [] = o_roots_Univariate_batch


arrExampleNumber = {'1','2','3','4'};
arrEmin = { 1e-12, 1e-11, 1e-10};
arrMeanMethod = {'None', 'Geometric Mean Matlab Method'};
arrEmax = 1e-12;
arrBoolAlphaTheta = {true, false};
arrLowRankApproxMethod = {'None', 'Standard STLN', 'Standard SNTLN'};
arrAPFMethod = {'None'};
arrDeconvolutionMethod = {'Separate', 'Batch', 'Batch Constrained'};
arrRoots_hx = { 'From ux', 'From Deconvolutions'};



global SETTINGS

SETTINGS.GCD_COEFFICIENT_METHOD = 'ux and vx';

for i0 = 1:1:length(arrRoots_hx)
    SETTINGS.ROOTS_HX = arrRoots_hx{i0};
    
    parfor i1 = 1:1:length(arrExampleNumber)
        ex_num = arrExampleNumber{i1};
        
        for i2 = 1:1:length(arrEmin)
            emin = arrEmin{i2};
            
            for i3 = 1:1:length(arrMeanMethod)
                mean_method = arrMeanMethod{i3};
                
                for i4 = 1:1:length(arrBoolAlphaTheta)
                    bool_alpha_theta = arrBoolAlphaTheta{i4};
                    
                    for i5 = 1:1:length(arrLowRankApproxMethod)
                        low_rank_approx_method = arrLowRankApproxMethod{i5};
                        
                        for i6 = 1:1:length(arrAPFMethod)
                            apf_method = arrAPFMethod{i6};
                            
                            for i7 = 1:1:length(arrDeconvolutionMethod)
                                deconvolution_method = arrDeconvolutionMethod{i7};
                                
                                file_name = 'log_roots.txt';
                                
                                try
                                    close all;
                                    clc;
                                    o_roots_Univariate(ex_num, emin, arrEmax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, deconvolution_method)
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






















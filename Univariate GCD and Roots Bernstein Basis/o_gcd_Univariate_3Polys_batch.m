function [] = o_gcd_Univariate_3Polys_batch()
% Perform a batch of three-polynomial GCD computations
%
% >> o_gcd_3Polys_batch
%
%


% Initialise array of example numbers
ex_num_arr = {'1','2','3','4','5','6','7','8','9','10','11','12'};

% Arrays of various noise levels
emin_arr = {1e-08,1e-10,1e-12};


% Initialise arrays related to preprocessing
bool_alpha_theta_arr = {true, false};
mean_method_arr = {'Geometric Mean Matlab Method', 'None'};


% Array of methods used to compute low rank approximation of the tth
% subresultant matrix - note that 'Standard STLN' and 'Standard SNTLN' are
% not working.
low_rank_approx_method_arr = {'None'};

% Array of methods for APF
apf_method_arr = {'None'};


% Only 
bool_log_arr = {false};

%
gcd_coefficient_method_arr = {'ux and vx'};

% Array of Sylvester subresultant matrix variants
%   'T'
%   'DT'
%   'TQ'
%   'DTQ'
Sylvester_Matrix_Variant_arr = {'T','DT','DTQ','TQ'};

% Array of Rank Revealing metrics
%   'Minimum Singular Values'
%
arrRankRevealingMetric = {'Minimum Singular Values'};

global SETTINGS



for i1 = 1:1:length(bool_log_arr)
    
    SETTINGS.BOOL_LOG = bool_log_arr{i1};
    
    for i2 = 1:1:length(gcd_coefficient_method_arr)
        
        SETTINGS.GCD_COEFFICIENT_METHOD = gcd_coefficient_method_arr{i2};
        
        % Changing example number
        parfor i3 = 1:1:length(ex_num_arr)
            
            ex_num = ex_num_arr{i3};
            
            % Changing lower noise boundary
            for i4 = 1:1:length(emin_arr)
                
                emin = emin_arr{i4};
                emax = 1e-12;
                
                % Changing low rank approx method
                for i5 = 1:1:length(low_rank_approx_method_arr)
                    
                    low_rank_approx_method = low_rank_approx_method_arr{i5};
                    
                    % Changing alpha theta boolean
                    for i6 = 1:1:length(bool_alpha_theta_arr)
                        
                        bool_alpha_theta = bool_alpha_theta_arr{i6};
                        
                        for i7 = 1:1:length(mean_method_arr)
                            mean_method = mean_method_arr{i7};
                            
                            for i8 = 1:1:length(Sylvester_Matrix_Variant_arr)
                                
                                sylvester_matrix_variant = Sylvester_Matrix_Variant_arr{i8};
                                
                                for i9 = 1:1:length(apf_method_arr)
                                    
                                    apf_method = apf_method_arr{i9};
                                    
                                    for i10 = 1: 1: length(arrRankRevealingMetric)
                                        
                                        rank_revealing_metric = arrRankRevealingMetric{i10};
                                        
                                        %try
                                            
                                            close all;
                                            clc;
                                            o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_matrix_variant, rank_revealing_metric)
                                            fileId = fopen('log_GCD_3Polys.txt','a')
                                            fprintf(fileId,'%s \n',datetime('now') , 'Success');
                                            fclose(fileId);
                                            
                                        %catch err
                                            
                                        %    fileId = fopen('log_GCD_3Polys.txt','a')
                                        %    fprintf(fileId,'%s %s \n\n\n',datetime('now'), getReport(err));
                                        %    fclose(fileId);
                                        %end
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




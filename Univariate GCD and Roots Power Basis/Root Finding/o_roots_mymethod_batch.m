function [] = o_roots_batch()

% Produce an array of Example numbers
%ex_num_arr = {'3'};

% Generate a set of examples
count = 1;
for i = 3:1:10
    ex_num_arr{count} = sprintf('Custom:m=%i low=-10 high=10',i);
    count = count + 1;
    
end

% Produce an array of lower noise level
emin_arr = ...
    {
    1e-12,...
    1e-11,...
    1e-10,...
    1e-9,...
    };

emax = 1e-12;

% Produce an array of mean methods
mean_method_arr = ...
    {...
    'None',...
    'Geometric Mean Matlab Method'
    };

% Produce an array to determine whether alpha and theta are computed.
bool_alpha_theta_arr = ...
    {
    'y',...
    'n'...
    };

%
low_rank_approx_method_arr = ...
    {
    'None',...
    'Standard STLN'
    };

%
deconvolution_method_arr = ...
    {
    'Batch',...
    'Separate'
    };

%
roots_ux_arr = {'From Deconvolution','From ux'};


global SETTINGS

for i6 = 1:1:length(roots_ux_arr)
    
    SETTINGS.ROOTS_UX = roots_ux_arr{i6};
    
    for i7 = 1:1:length(deconvolution_method_arr)
        
        SETTINGS.DECONVOLUTION_METHOD = deconvolution_method_arr{i7};
        
        parfor i1 = 1:1:length(ex_num_arr)
            
            ex_num = ex_num_arr{i1};
            
            for i2 = 1:1:length(emin_arr)
                
                emin = emin_arr{i2};
                
                for i3 = 1:1:length(mean_method_arr)
                    
                    mean_method = mean_method_arr{i3};
                    
                    for i4 = 1:1:length(bool_alpha_theta_arr)
                        
                        bool_alpha_theta = bool_alpha_theta_arr{i4};
                        
                        for i5 = 1:1:length(low_rank_approx_method_arr)
                            
                            low_rank_approx_method = low_rank_approx_method_arr{i5};

                            try
                                
                                o_roots(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method);
                            catch err
                                fprintf(err.message);
                            end
                        end
                    end
                end
            end
        end
    end
    
end

end
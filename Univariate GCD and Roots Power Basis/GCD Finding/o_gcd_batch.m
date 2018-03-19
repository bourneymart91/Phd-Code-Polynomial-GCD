function []  = o_gcd_batch()

%
%ex_num_arr = {'1','2','3','4','5','6','7','8','9','10','11','12'};
% Custom:m=(\d+).n=(\d+).t=(\d).low=(-?\d+).high=(-?\d+)

count = 1;
for i = 3:1:10
    for j = 3:1:i
        for k = 1:1:j
            ex_num_arr{count} = sprintf('Custom:m=%i n=%i t=%i low=-1 high=5',i,j,k);
            count = count +1;
        end
    end
end

% Set upper noise to be 1e-8
el_arr = ...
    {
     %1e-8,...
     %1e-9,...
     %1e-10,...
     %1e-11,...
     1e-12
    };

emin = 1e-12;


% mean method
mean_method_arr = ...
    {
    'None',...
    'Geometric Mean Matlab Method'
    };

%
alpha_theta_arr = ...
    {
    'y',...
    'n'
    };

%
low_rank_approx_method_arr = ...
    {
    'None',...
    'Standard STLN'
    };

vars = [];
count = 1;
for i1 = 1:1:length(ex_num_arr)
    for i2 = 1:1:length(mean_method_arr)
        for i3 = 1:1:length(el_arr)
            for i4 = 1:1:length(alpha_theta_arr)
                for i5 = 1:1:length(low_rank_approx_method_arr)                    
                    vars(count,:) = [i1,i2,i3,i4,i5];
                    count = count + 1;
                end
            end
        end
    end
end

emax = 1e-12;
parfor i = 1:1:size(vars,1)
    
    var = vars(i,:);
    ex_num = ex_num_arr{var(1)};
    mean_method = mean_method_arr{var(2)};
    emin = el_arr{var(3)};
    bool_alpha_theta = alpha_theta_arr{var(4)};
    low_rank_approx_method = low_rank_approx_method_arr{var(5)};
    
    try
        o_gcd_2Polys(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method);
    catch err
        fprintf(err.message);
    end
    
end


end
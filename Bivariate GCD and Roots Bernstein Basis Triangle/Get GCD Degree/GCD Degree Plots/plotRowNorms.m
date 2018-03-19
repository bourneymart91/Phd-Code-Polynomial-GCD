function [] = plotRowNorms(arr_RowNorms, limits_k, limits_t)
%
% % Inputs
%
% arr_RowNorms : (Array)
%
% limits_k : (Int Int)
%
% limits_t : (Int Int)

% Get limits
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);



figure_name = 'Plotting Row Norms';
figure('name',figure_name)
hold on

for i = lowerLimit_k:1:upperLimit_k
    
    vRowNorms = arr_RowNorms{i};
    
    vec_i = i.*ones(length(vRowNorms),1);
   
    plot(vec_i, vRowNorms); 
    
end

% Plot vertical lines
vline(limits_t, {'-r','-r'});


hold off

end
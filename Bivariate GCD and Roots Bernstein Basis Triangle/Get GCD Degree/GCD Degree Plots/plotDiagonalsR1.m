function [] = plotDiagonalsR1(arr_R1, limits_k, limits_t)
%
% % Inputs
%
% arr_R1 : (Array)
%
% limits_k : (Int Int)
%
% limits_t : (Int Int)
%
% 

% Get my limits
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

% Get actual limits
lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

% Get number of subresultants constructed
nSubresultants = upperLimit_k - lowerLimit_k + 1;

global SETTINGS

figure_name = sprintf('Diagonals of R1 from QR decomposition of %s',SETTINGS.SYLVESTER_MATRIX_VARIANT);
figure('name',figure_name)
hold on

for i = 1 : 1 : nSubresultants
   
    k = lowerLimit_k + (i-1);
    
    % Get diagonal entries of R1
    temp_vec = diag(arr_R1{i});
    
    % Get vector of [ i i i i ...]
    vec_x = k*ones(length(temp_vec));
    
    % plot
    plot(vec_x, log10(temp_vec), '*')
    
    
end

% Add vertical lines to show limits
vline(lowerLimit_t);
vline(upperLimit_t);
hold off

end
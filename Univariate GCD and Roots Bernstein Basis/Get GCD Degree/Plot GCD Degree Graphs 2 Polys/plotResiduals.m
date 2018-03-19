function plotResiduals(vMinimumResiduals, limits_k, limits_t, rank_range)
% Plot the minimum residuals associated with removal of the optimal column
% from each subresultant matrix S_{k}(f,g)
%
% % Inputs
%
% vMinimumResiduals : (Vector) Vector of residuals
%
% limits_k : [Int Int] Limits on k
%
% limits_t : [Int Int] Limits on possible t values
%
% rank_range : [Float, Float] 


lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

% % Plot Residuals
global SETTINGS

figure_name = sprintf('%s : Residuals of %s', mfilename, ...
    SETTINGS.SYLVESTER_MATRIX_VARIANT);

figure('name', figure_name);

hold on
vec_x = lowerLimit_k:1:upperLimit_k;
plot(vec_x, log10(vMinimumResiduals), '-s', 'DisplayName', 'Residuals by SVD')
ylabel('log_{10} Residuals')
xlabel('k')
legend(gca,'show');

hline([rank_range mean(rank_range)],{'-r','-r','-b'})
vline(limits_t,{'r','r'})


hold off


end






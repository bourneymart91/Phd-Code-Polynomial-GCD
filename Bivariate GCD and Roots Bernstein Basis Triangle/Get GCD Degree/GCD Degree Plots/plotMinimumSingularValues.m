function [] = plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t, rank_range)
%
% % Inputs
%
% vMinimumSingularValues : (Vector)
%
% limits_k : (Int Int)
%
% limits_t : (Int Int)
%
% rank_range : [Float Float]

%
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

vec_x = lowerLimit_k : 1 : upperLimit_k;
global SETTINGS

figure_name = sprintf('Minimum Singular Values of %s' ,SETTINGS.SYLVESTER_MATRIX_VARIANT);
figure('name',figure_name)
hold on

plot(vec_x, log10(vMinimumSingularValues),'-s','LineWidth',2)

try
    xlim(limits_k);
catch
end

% 
hline(rank_range,{'--r','--r'});
vline(limits_t, {'--r','--r'});


% Labels
xlabel('$k$', 'Interpreter', 'latex')
ylabel('$\log_{10}\left( \sigma_{k} \right)$', 'Interpreter', 'latex')

% Display
grid on
box on





hold off

end